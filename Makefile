SHELL := bash -O extglob
all:
	cc -std=c11 -pedantic -g -Wall -lm -D_XOPEN_SOURCE=600 -o WeberSim main.c
clean:
	rm !(main.c||Makefile||.git||.||..)
