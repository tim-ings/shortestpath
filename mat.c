// Tim Ings 21716194
#include "mat.h"


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
int* matNew(int n) {
    return mallocRetry(n * n * sizeof(int));
}

// Gets a value from matrix mat at (i,j) with dimensions n
int matGet(int* mat, int n, int i, int j) {
    return mat[i * n + j];
}

// Sets the value at (i,j) to v in matrix mat with dimensions n
void matSet(int* mat, int n, int i, int j, int v) {
    mat[i * n + j] = v;
}

// Prints a matrix to stdout
// Does not line up cols when integers get big
void matPrint(int *mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", matGet(mat, n, i, j));
        }
        printf("\n");
    }
}

// Loads an adjacency matrix from binary file
// IN: filepath
// OUT: n
int* matFromBinFile(const char* filepath, int* n) {
    FILE* fp;
    fp = fopen(filepath, "rb");
    fread(n, sizeof(int), 1, fp);
    int* mat = matNew(*n);
    fread(mat, sizeof(int), (*n) * (*n), fp);
    // convert all zeroes that are not i==j in input mat to INT_INF
    for (int i = 0; i < (*n); i++) {
        for (int j = 0; j < (*n); j++) {
            int val = matGet(mat, (*n), i, j);
            if (val == 0 && i != j) {
                matSet(mat, (*n), i, j, INT_INF);
            } else if (i == j) {
                matSet(mat, (*n), i, j, 0);
            }
        }
    }
    fclose(fp);
    return mat;
}

// this assumes that INT_INF has been converted to 0 already
void matToBinFile(const char* filepath, int* mat, int n) {
    FILE* fp;
    fp = fopen(filepath, "wb");
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(mat, sizeof(int), n * n, fp);
    fclose(fp);
}
