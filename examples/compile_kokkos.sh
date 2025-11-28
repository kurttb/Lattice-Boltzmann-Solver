#!/bin/bash
KOKKOS_PATH=/Users/ethanstout/packages/kokkos
clang++ -I${KOKKOS_PATH}/include -I../include -L${KOKKOS_PATH}/lib -lkokkoscore -fopenmp -O3 -march=native -ffast-math naive_kokkos.cpp ../src/VtkWriter.cpp -o prog

export OMP_PROC_BIND=spread  
export OMP_PLACES=threads
