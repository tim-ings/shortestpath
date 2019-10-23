#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "mat.h"

#define RANK_ROOT 0

const int INFINITY = INT_MAX / 2;

void floydWarshall(int* partMat, int nMain, int commRank, int commSize, MPI_Comm comm) {
    int root;
    int* kRow = malloc(nMain * sizeof(int));
    int kMain;
    for (kMain = 0; kMain < nMain; kMain++) {
        int kRank = kMain / (nMain / commSize);
        // copy row k into kRow
        if (commRank == kRank) {
            int kPart = kMain % (nMain / commSize);
            int kRowJ;
            for (kRowJ = 0; kRowJ < nMain; kRowJ++) {
                kRow[kRowJ] = matGet(partMat, nMain, kPart, kRowJ);
            }
        }
        MPI_Bcast(kRow, nMain, MPI_INT, kRank, comm);
        int iPart;
        for (iPart = 0; iPart < nMain / commSize; iPart++) {
            int jMain;
            for (jMain = 0; jMain < nMain; jMain++) {
                int val = matGet(partMat, nMain, iPart, kMain);
                val += kRow[jMain];
                if (val < matGet(partMat, nMain, iPart, jMain)) {
                    matSet(partMat, nMain, iPart, jMain, val);
                }
            }
        }
    }
    free(kRow);
}

int main(int argc, char* argv[]) {
    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);

    // Returns the size of the group associated with a communicator
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // Determines the rank of the calling process in the communicator
    int commRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    if (argc < 2 && commRank == RANK_ROOT) {
        printf("Usage: %s <file> [-o] [-t] [-p]\n"
               "\t<file>\tBinary file containg an adjacency matrix\n"
               "\t-o\tIndicates that an output file should be written\n"
               "\t-t\tIndicates that the operation should be timed\n"
               "\t-p\tIndicates that the result should be printed\n", argv[0]);
        // Terminates MPI execution environment
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit(EXIT_FAILURE);
    }

    // Load adjacency matrix from binary file
    int inN;
    int* inMat;
    if (commRank == RANK_ROOT) {
        inMat = matFromBinFile(argv[1], &inN);
    }

    // Broadcasts a message from the process with rank root to all other processes of the group
    MPI_Bcast(&inN, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Sends data from one task to all tasks in a group
    int partN = inN * inN / commSize;
    int* partMat = matNew(partN);
    MPI_Scatter(inMat, partN, MPI_INT, partMat, partN, MPI_INT, RANK_ROOT, MPI_COMM_WORLD);

    double t0 = MPI_Wtime();
    floydWarshall(partMat, inN, commRank, commSize, MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    // Gathers values from a group of processes
    int finalN = inN;
    int* finalMat = matNew(finalN);
    MPI_Gather(partMat, partN, MPI_INT, finalMat, partN, MPI_INT, RANK_ROOT, MPI_COMM_WORLD);

    if (commRank == RANK_ROOT) {
        matPrint(finalMat, finalN);
        char outpath[BUFSIZ];
        strcat(outpath, argv[1]);
        strcat(outpath, ".out");
        matToBinFile(outpath, finalMat, finalN);
    }

    // Terminates MPI execution environment
    free(partMat);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
