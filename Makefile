SHELL := bash -O extglob
all:
	cc -Ofast -fopenmp -std=c18 -pedantic -Wall -D_XOPEN_SOURCE=600 -o WeberSim main.c -lm
debug:
	cc -g -pg -Ofast -fopenmp -std=c18 -pedantic -Wall -D_XOPEN_SOURCE=600 -o WeberSim main.c -lm
clean:
	rm !(*.c||*.h||Makefile||*.blend||.git||.||..||*.py||WeberSim||*.gnuplot||*.tex)
