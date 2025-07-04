#!/bin/bash

## run this to compile without optimisation.
gcc -I/home/lmj/R/R-4.4.1/include -DNDEBUG -I/usr/local/include -fpic -g -Og -c read_bam.c -o read_bam.o
gcc -shared -L/usr/local/lib -o read_bam.so read_bam.o common.o -lhts
