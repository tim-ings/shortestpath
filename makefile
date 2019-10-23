default: build

build:
	mpicc -std=c99 -lm *.c -o main.out
