// Tim Ings 21716194
// This file defines functions that are used to manipulate matrices
#include "mat.h"


// wrapper around malloc that allows us to retry a failed malloc up to MALLOC_MAX_RETRY times
void* mallocRetry(size_t s) {
#ifdef MALLOC_SHOULD_RETRY
    void* mem = NULL;
    int attempts = 0;
    while (mem == NULL && attempts < MALLOC_MAX_RETRY) {
        mem = malloc(s);
        attempts++;
    }
    if (mem == NULL) {
        printf("WARNING: FAILED TO ALLOCATE MEMORY AFTER %d RETRIES\n", attempts);
        exit(EXIT_FAILURE);
    }
    return mem;
#else
    return malloc(s);
#endif
}

// Allocates a new matrix
// we add 1 to the requested size to avoid segfaults with larger matrices (over ~512)
int* matNew(int n) {
    int nn = n + 1;
    return mallocRetry(nn * nn * sizeof(int));
}

// Gets a value from matrix mat at (i,j) with size n
int matGet(int* mat, int n, int i, int j) {
    return mat[i * n + j];
}

// Sets the value at (i,j) to v in matrix mat with size n
void matSet(int* mat, int n, int i, int j, int v) {
    mat[i * n + j] = v;
}

// Prints a matrix to stdout
// Does not line up cols when integers get big, mostly used for debugging
// could be used to output text files
void matPrint(int *mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matGet(mat, n, i, j));
        }
        printf("\n");
    }
}

// Loads an adjacency matrix from binary file
// we output the size of the matrix to n
int* matFromBinFile(const char* filepath, int* n) {
    FILE* fp;
    fp = fopen(filepath, "rb");
    fread(n, sizeof(int), 1, fp); // read the mat size
    int* mat = matNew(*n);
    fread(mat, sizeof(int), (*n) * (*n), fp); // read the matrix values for all (i,j)
    
    for (int i = 0; i < (*n); i++) {
        for (int j = 0; j < (*n); j++) {
            int val = matGet(mat, (*n), i, j);
            if (val == 0 && i != j) {
                // convert all zeroes that are not i==j in input mat to INT_INF
                matSet(mat, (*n), i, j, INT_INF);
            } else if (i == j) {
                // ensure the diagonal is all 0
                matSet(mat, (*n), i, j, 0);
            }
        }
    }
    fclose(fp);
    return mat;
}

void matToBinFile(const char* filepath, int* mat, int n) {
    // convert INT_INF to zero for output
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int val = matGet(mat, n, i, j);
            if (val == INT_INF) {
                matSet(mat, n, i, j, 0);
            }
        }
    }
    FILE* fp;
    fp = fopen(filepath, "wb");
    fwrite(&n, sizeof(int), 1, fp); // write the size of the matrix
    fwrite(mat, sizeof(int), n * n, fp); // write values for all (i,j)
    fclose(fp);
}
