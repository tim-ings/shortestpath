// Tim Ings 21716194
// This file declares functions that are used to manipulate matrices
//
// if MALLOC_SHOULD_RETRY is defined then all calls to malloc 
// will be repeated MALLOC_MAX_RETRY times. This is useful for larger
// matrices as sometimes the clusters do not have enough memory
#ifndef MAT_H
#define MAT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define INT_INF INT_MAX / 2
#define MALLOC_MAX_RETRY 50
#define MALLOC_SHOULD_RETRY


void* mallocRetry(size_t s);
int* matNew(int n);
int matGet(int* mat, int n, int i, int j);
void matSet(int* mat, int n, int i, int j, int v);
void matPrint(int *mat, int n);
int* matFromBinFile(const char* filepath, int* n);
void matToBinFile(const char* filepath, int* mat, int n);

#endif