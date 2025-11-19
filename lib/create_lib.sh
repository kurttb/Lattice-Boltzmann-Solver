#!/bin/bash

# Compile each src file
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/BoundaryConditions.cpp -c
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/Collisions.cpp -c
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/ComputeState.cpp -c
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/D2Q9Problem.cpp -c
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/Streaming.cpp -c
clang++ -I../include -O3 -ffast-math -march=native -fopenmp ../src/VtkWriter.cpp -c


clang++ -dynamiclib -fopenmp -O3 -ffast-math -march=native -o libLBM.dylib *.o # MacOS
#g++ -shared -fopenmp -O3 -ffast-math -march=native -o libLBM.so *.o #Linux
#ar rcs libLBM.a *.o
rm *.o
