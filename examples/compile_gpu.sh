#!/bin/bash 
KOKKOS_PATH=$HOME/.local/kokkos 
CXX=${KOKKOS_PATH}/bin/nvcc_wrapper 

$CXX naive_kokkos.cpp ../src/VtkWriter.cpp --expt-extended-lambda -std=c++17 -I${KOKKOS_PATH}/include -I../include -L${KOKKOS_PATH}/lib64 -lkokkoscore -lcuda -fopenmp -O3 -march=native -ffast-math -arch=sm_70 -o prog
