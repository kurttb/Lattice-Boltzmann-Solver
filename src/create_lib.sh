#!/bin/bash
g++ -I../include -c vtk_writer.cpp -o vtk_writer.o
ar rcs libLBM.a vtk_writer.o
#g++ -dynamiclib -o libLBM.dylib vtk_writer.o
rm vtk_writer.o
