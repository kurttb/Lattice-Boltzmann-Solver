#!/bin/bash
KOKKOS_PATH=${HOME}/packages/kokkos
clang++ -std=c++17 -I${KOKKOS_PATH}/include -I../include -L${KOKKOS_PATH}/lib Couette_kokkos_case.cpp ../src/D2Q9Problem.cpp -lkokkoscore -fopenmp -O3 -march=native -ffast-math -o prog

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
