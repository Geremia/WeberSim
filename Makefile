SHELL := bash -O extglob
all:
	cc -Ofast -std=c11 -pedantic -Wall -lm -D_XOPEN_SOURCE=600 -o WeberSim main.c
debug:
	cc -g -Ofast -std=c11 -pedantic -Wall -lm -D_XOPEN_SOURCE=600 -o WeberSim main.c
clean:
	rm !(*.c||*.h||Makefile||.git||.||..||*.py)
