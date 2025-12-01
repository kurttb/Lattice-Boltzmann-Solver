#!/bin/bash 
KOKKOS_PATH=$HOME/.local/kokkos
CXX=${KOKKOS_PATH}/bin/nvcc_wrapper

$CXX Couette_kokkos_case.cpp ../src/D2Q9Problem.cpp ../src/VtkWriter.cpp --expt-extended-lambda -std=c++17 -I${KOKKOS_PATH}/include -I../include -L${KOKKOS_PATH}/lib64 -lkokkoscore -lcuda -O3 -ffast-math -arch=sm_70 -o prog
