/* 
Name:
    Tim Ings 21716194

Compile:
    make
  -OR-
    mpicc -std=c99 *.c -o main.out

Run:
    mpirun [-np <numprocs>] [--hostfile <hostfile>] main.out <file> [-p] [-o] [-t]
    
    Required:
        <file>  path to binary file containing input adjacency matrix
    Optional:
        --hostfile  Specify MPI hostfile
        -np         Number is processes to launch
        -p          Prints final adjacency matrix to stdout 
        -o          Outputs matrix to binary file
        -t          Times the operation and outputs result to stdout

The number of vertices in the input matrix should be evenly divisible by the number of processes requested

Compiling with DEBUG defined will print out the next action that program will take on the root process
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "mat.h"

#define RANK_ROOT 0

// simple bool implementation for pre c99 compatibility
#define bool int
#define true 1
#define false 0

// Compiling with DEBUG defined will print out the next action that program will take on the root process
//#define DEBUG

// struct to hold parsed command line args
struct {
    bool shouldOutput;
    bool shouldTime;
    bool shouldPrint;
    char* filepath;
} typedef CLARGS;


// floyd warshall implementation that works on part mats and utilises Open MPI
void floydWarshall(int* partMat, int nMain, int commRank, int commSize, MPI_Comm comm) {
    int root;
    // malloc an array to store row k
    int* kRow = mallocRetry(nMain * sizeof(int));
    int kMain; // pre c99 compatibility
    for (kMain = 0; kMain < nMain; kMain++) {
        int kRank = kMain / (nMain / commSize); // get the rank of the owner of this row
        // copy row k into kRow if this is our processes' k row
        if (commRank == kRank) {
            int kPart = kMain % (nMain / commSize);
            int kRowJ;
            for (kRowJ = 0; kRowJ < nMain; kRowJ++) {
                kRow[kRowJ] = matGet(partMat, nMain, kPart, kRowJ);
            }
        }
        // Broadcasts a message from the process with rank root to all other processes of the group
        // We boadcast the row k to all other processes using kRank as root
        // that is, the rank of the process who has been assinged this row by MPI_Scatter
        MPI_Bcast(kRow, nMain, MPI_INT, kRank, comm);
        int iPart;
        // update shorest paths
        for (iPart = 0; iPart < nMain / commSize; iPart++) {
            int jMain;
            for (jMain = 0; jMain < nMain; jMain++) {
                // calculate a potentially shorter path
                int val = matGet(partMat, nMain, iPart, kMain);
                val += kRow[jMain];
                // if it is shorter, then update our matrix
                if (val < matGet(partMat, nMain, iPart, jMain)) {
                    matSet(partMat, nMain, iPart, jMain, val);
                }
            }
        }
    }
    free(kRow);
}

// prints applications usage to stdout
void printUsage(const char* name) {
    printf("Usage: %s <file> [-o] [-t] [-p]\n"
           "\t<file>\tPath to binary file containg an adjacency matrix\n"
           "\t-o\tOutputs to binary file in-path.out\n"
           "\t-t\tTimes the operation\n"
           "\t-p\tPrints result to stdout\n", name);
}

// parses command line arguments
// returns false on invalid argc/argv
// outputs to argsOut
bool parseArgs(int argc, char* argv[], CLARGS* argsOut) {
    // get the file path or fail
    if (argc < 2) {
        return false;
    }
    argsOut->filepath = argv[1];
    // get every optional arg
    for (int i = 2; i < argc; i++) {
        char* v = argv[i];
        if (strcmp(v, "-o") == 0) {
            argsOut->shouldOutput = true;
        } else if (strcmp(v, "-t") == 0) {
            argsOut->shouldTime = true;
        } else if (strcmp(v, "-p") == 0) {
            argsOut->shouldPrint = true;
        }
    }
    return true;
}

// entry point
int main(int argc, char* argv[]) {
    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);

    // Returns the size of the group associated with a communicator
    // we use this to partition our input matrix
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // Determines the rank of the calling process in the communicator
    // we use this to identify the current process
    int commRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START PARSING CLI ARGS\n");
#endif
    // parse command line args into args structure
    CLARGS args;
    if (!parseArgs(argc, argv, &args) && commRank == RANK_ROOT) {
        printUsage(argv[0]);
        // Terminates MPI execution environment
        // we do this to signal to all out non root processes to terminate
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START FILE LOADING ON ROOT\n");
#endif
    // Load adjacency matrix from binary file
    int inN;
    int* inMat;
    if (commRank == RANK_ROOT) {
        // load the matrix from binary file on the root process
        inMat = matFromBinFile(args.filepath, &inN);
        if (inN % commSize != 0) {
            // TODO: could adjust the MPI_Group size here to next highest valid commSize instead of aborting
            // Would need to keep track of our new comm as a variable instead of using MPI_COMM_WORLD
            // MPI_Comm_group -> MPI_Group_range_excl -> MPI_Comm_create
            // We require this as we do not want to split up rows over multiple processes
            // See MPI_Bcast() call in floydWarshall()
            printf("The number of vertices in the input matrix "
                   "must be evenly divisible by the number of "
                   "processes requested.\n");
            // Terminates MPI execution environment
            // we do this to signal to all out non root processes to terminate
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            exit(EXIT_FAILURE);
        }
    }

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO BCASR IN N\n");
#endif
    // Broadcasts a message from the process with rank root to all other processes of the group
    // broadcast size of input matrix to all non root processes
    MPI_Bcast(&inN, 1, MPI_INT, RANK_ROOT, MPI_COMM_WORLD);

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START TIMING\n");
#endif
    // Synchronization between MPI processes in a group
    // we do this to ensure all procs are up to this point before we start timing
    MPI_Barrier(MPI_COMM_WORLD); 
    // Returns an elapsed time on the calling processor
    // get the start time of our algorithm
    double t0 = MPI_Wtime();

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START SCATTERING\n");
#endif

    // Sends data from one task to all tasks in a group
    // partition the input matrix by scattering the it to all our processes
    int partN = inN * inN / commSize;
    int* partMat = matNew(partN);
    MPI_Scatter(inMat, partN, MPI_INT, partMat, partN, MPI_INT, RANK_ROOT, MPI_COMM_WORLD);

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START FLOYD\n");
#endif
    // run the floyd-warshall algorithm on our part matrix
    floydWarshall(partMat, inN, commRank, commSize, MPI_COMM_WORLD);
    
#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO START GATHERING\n");
#endif
    // Gathers values from a group of processes
    // this is the opposite of scatter
    // we use this to reconstruct our final matrix
    int finalN = inN;
    int* finalMat = matNew(finalN);
    MPI_Gather(partMat, partN, MPI_INT, finalMat, partN, MPI_INT, RANK_ROOT, MPI_COMM_WORLD);

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO STOP TIMING\n");
#endif
    // Returns an elapsed time on the calling processor
    // get the end time of our algorithm
    double t1 = MPI_Wtime();

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO OUTPUT REQUIRED DATA\n");
#endif
    // Outputs the requested data on root
    if (commRank == RANK_ROOT) {
        if (args.shouldPrint) {
            matPrint(finalMat, finalN);
        }
        if (args.shouldTime) {
            printf("Floyd-Warshall running on %d processes for matrix size %d took %.4fms\n", commSize, inN, (t1 - t0) * 1000);
        }
        if (args.shouldOutput) {
            char outpath[BUFSIZ];
            strcat(outpath, argv[1]);
            strcat(outpath, ".out");
            matToBinFile(outpath, finalMat, finalN);
        }
    }

#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO FREE MEMORY\n");
#endif
    // free malloc'd memory
    free(inMat);
    free(partMat);
    free(finalMat);
#ifdef DEBUG
    if (commRank == RANK_ROOT) printf("\tABOUT TO FINALISE\n");
#endif
    // Terminates MPI execution environment
    // we use this to exit gracefully
    MPI_Finalize();
    return EXIT_SUCCESS;
}
