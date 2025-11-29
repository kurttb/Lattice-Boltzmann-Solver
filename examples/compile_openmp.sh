#!/bin/bash
KOKKOS_PATH=${HOME}/.local/kokkos
g++ -std=c++17 -I${KOKKOS_PATH}/include -I../include -L${KOKKOS_PATH}/lib naive_kokkos.cpp ../src/VtkWriter.cpp -lkokkoscore -ldl -fopenmp -O3 -march=native -ffast-math -o prog

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
