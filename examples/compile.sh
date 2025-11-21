#!/bin/bash

clang++ -O3 -ffast-math -march=native -fopenmp -I../include -L../lib -lLBM Couette.cpp -o prog
