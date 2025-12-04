# Lattice-Boltzmann-Solver

This project implements a modular, high-performance Lattice Boltzmann Method (LBM) solver for 2D flows on general rectangular domains.
Built with Kokkos, the library targets both CPU and GPU architectures from a single codebase, allowing for portable performance across modern HPC systems.

**Highlights**
- Single codebase runs efficiently on multi-core CPUs and NVIDIA GPUs
- Verified on Poiseuille channel, planar Couette, and lid-driven cavity flows
- Arbitrary combinations of bounce-back, velocity, and periodic boundary conditions
- Simple 15–20 line user code for new cases

### Performance Summary (single-precision, full simulation with boundary conditions)

| Hardware                                   | Threads / Device | MLUPS   | Notes                                      |
|--------------------------------------------|------------------|---------|--------------------------------------------|
| Great Lakes (Intel Xeon Gold 6154, Skylake) | 32 threads       | **527** | Bandwidth-saturated on a single socket     |
| Tesla V100 (NERSC Perlmutter / Great Lakes GPU) | 1 GPU            | **3520**| Kokkos CUDA                                |
| RTX 5070 Ti (personal workstation)        | 1 GPU            | **5500**| Kokkos CUDA                                |

## Downloading the solver
The solver can be downloaded by cloning it using the command:
```
git clone https://github.com/kurttb/Lattice-Boltzmann-Solver.git
```

## Building the Software

The instructions for building the solver are as follows:

### Dependencies
The solver requires several dependencies in order to build the software. Different dependencies are required depending on the architecture in which you intend to build for. For all builds, cmake and Kokkos are required. Paraview is also required for visualization. A deprecated pure OpenMP branch exists (modular_omp), but is no longer maintained
```
git clone -b modular_omp https://github.com/kurttb/Lattice-Boltzmann-Solver.git 
```
For Kokkos builds, the required backends must be built on your machine. Using an OpenMP backend requires having OpenMP built on your machine, and using a CUDA backend requires the NVIDIA CUDA toolkit. Instructions for building Kokkos can be found at:
```
https://github.com/kokkos
```
When building Kokkos, make sure to specify the build backend and the CPU/GPU architecture you wish to run on to get maximum performance optimizations

### Configuring, Compiling and Installing
The build can be configured using cmake and by passing in the appropriate options:
```
cd $PROJECT_ROOT
mkdir build && cd build
cmake -OPTIONS ../
```
Common configuration options include:
	-DKokkos_ROOT			- Specifies the root directory in which Kokkos is installed. 
							  By default it is set to /usr/local, but if Kokkos is installed 
							  elsewhere you must specify the path to the instillaiton.
	-DCMAKE_CXX_COMPILER	- Specifies the c++ compiler to be used in the build. For builds with
							  the CUDA backend, use ${DKokkos\_ROOT}/bin/nvcc\_wrapper}
	-DCMAKE_INSTALL_PREFIX  - Set the installation prefix.
	-SHARED_LIB				- Set to ON if you wish to build a shared library instead of a static
							  library (a static library will be built by default)
	-DCMAKE_BUILD_TYPE		- By default, the build type is set to release. However, debug builds 
							  can be set by passing in Debug to this configuration option

After configuring, the project can be compiled and installed in the build directory using:
```
make
make install
```
The core library will by default be built in ${PROJECT\_ROOT}/lib, and the example exectuables
will be placed in ${PROJECT\_ROOT}/bin. You can add a case of your own by naming the case
${CASE\_NAME}\_case.cpp and by placing this file in ${PROJECT\_ROOT}/examples.

### Additional Useful Commands
	make clean
		Removes the library and the built exectuables but keeps the configuration options.
	
	make uninstall
		Removes all installed files from 'make install'.

## Examples
Several examples exist in ${PROJECT\_ROOT}/examples. Cases include Couette Flow, a Lid-Driven cavity, and Poiseuille Flow.

