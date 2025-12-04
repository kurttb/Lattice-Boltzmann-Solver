# Lattice-Boltzmann-Solver

This project implements a modular, high-performance Lattice Boltzmann Method (LBM) solver for 2D flows on general rectangular domains.
Built with Kokkos, the library targets both CPU and GPU architectures from a single codebase, allowing for portable performance across modern HPC systems.

**Highlights**
- Single codebase runs efficiently on multi-core CPUs and NVIDIA GPUs
- Modular, object-oriented design
- Easily extendable to incorporate additional physics (i.e. forcing, turbulence models)
- Verified on Poiseuille channel, planar Couette, and lid-driven cavity flows
- Arbitrary combinations of bounce-back, velocity, and periodic boundary conditions
- Simple 15–20 line user code for new cases

### Performance Summary (single-precision, full simulation with boundary conditions)

| Hardware                                   | Threads / Device | MLUPS   | Notes                                      |
|--------------------------------------------|------------------|---------|--------------------------------------------|
| Great Lakes (Intel Xeon Gold 6154, Skylake) | 32 threads       | **527** | Bandwidth-saturated on a single socket     |
| Tesla V100 (Michigan Great Lakes GPU) | 1 GPU            | **3520**| Kokkos CUDA                                |
| RTX 5070 Ti (personal workstation)        | 1 GPU            | **5500**| Kokkos CUDA                                |

## Downloading the solver
The solver can be downloaded by cloning it using the command:
```
git clone https://github.com/kurttb/Lattice-Boltzmann-Solver.git
```

A deprecated pure OpenMP branch exists (modular_omp), but is no longer maintained:
```
git clone -b modular_omp https://github.com/kurttb/Lattice-Boltzmann-Solver.git 
```

## Building the Software

The instructions for building the solver are as follows:

### Dependencies
Different dependencies are required depending on the backend (Serial/OpenMP/CUDA) that you intend to build. For all builds, cmake and Kokkos are required. Paraview is also recommended for visualization. For Kokkos builds, the required backends must be built on your machine. Using an OpenMP backend requires having OpenMP built on your machine, and using a CUDA backend requires the NVIDIA CUDA toolkit. Instructions for building Kokkos can be found at:
```
https://github.com/kokkos
```
When building Kokkos, you must specify the build backend (Serial/OpenMP/CUDA), as well as the CPU/GPU architecture you wish to run on to obtain maximum performance optimizations

### Configuring, Compiling and Installing
The build can be configured using cmake and by passing in the appropriate options:
```
cd $PROJECT_ROOT
mkdir build && cd build
cmake -OPTIONS ../
```
### Common CMake Options

| Option                         | Description                                                                                      |
|---------------------------     |--------------------------------------------------------------------------------------------------|
| `-DKokkos_ROOT=/path`          | Path to Kokkos installation (default: `/usr/local`)                                              |
| `-DCMAKE_CXX_COMPILER`         | C++ compiler to use. For CUDA builds, set to `${Kokkos_ROOT}/bin/nvcc_wrapper`                   |
| `-DCMAKE_INSTALL_PREFIX`       | Installation directory (if you run `make install`)                                               |
| `-DBUILD_SHARED_LIBS=ON`       | Build shared library instead of static (default: static)                                         |
| `-DCMAKE_BUILD_TYPE`           | `Release` (default, fast) or `Debug`                                                             |

**Example (OpenMP build):**
```bash
cmake    -DKokkos_ROOT=$HOME/kokkos-openmp-install \
         -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
         -DCMAKE_BUILD_TYPE=Release \
         ../
```

**Example (CUDA build):**
```bash
cmake    -DKokkos_ROOT=$HOME/kokkos-cuda-install \
         -DCMAKE_CXX_COMPILER=$HOME/kokkos-install/bin/nvcc_wrapper \
         -DCMAKE_BUILD_TYPE=Release \
         ../
```

After configuring, the project can be compiled and installed in the build directory using:
```
make
make install
```
Binaries appear in `bin/`, the library in `lib/`.

Add your own simulation case by placing a file named `CASE_NAME_case.cpp` in the `examples/` folder. It will be automatically built.

### Additional Useful Commands
	make clean
		Removes the library and the built exectuables but keeps the configuration options.
	
	make uninstall
		Removes all installed files from 'make install'.

## Examples
Several examples exist in examples/. Cases include Couette Flow, a Lid-Driven cavity, and Poiseuille Flow. A minimal example for Couette flow is shown below:

**Couette Flow Example:**
``` c++
#include <string>
#include <cmath>
#include "D2Q9Problem.hpp"

int main() {

    // Grid resolution
    const int Nx = 100;
    const int Ny = 100;

    // Flow parameters
    const float Re = 100.0f;   // Reynolds number
    const float Ma = 0.1f;     // Mach number
    const int   LChar = Ny;    // Characteristic length scale

    // Derived quantities
    float cs    = 1.0f / std::sqrt(3.0f);
    float uChar = Ma * cs;
    float nu    = uChar * LChar / Re;   // Kinematic viscosity

    // Create LBM problem
    auto prob = LBM::D2Q9Problem(Nx, Ny);

    prob.setViscosity(nu);
    prob.setIC(1.0f, 0.0f, 0.0f);        // rho0, ux0, uy0
    prob.setNumTimeSteps(50000);

    // Boundary conditions
    prob.setBC("Top",    "WallTangentVelocity", uChar);
    prob.setBC("Bottom", "BounceBack");
    prob.setBC("Left",   "Periodic");
    prob.setBC("Right",  "Periodic");

    // Body force
    prob.setForces(0.0f, 0.0f);

    // Run and write output
    prob.runSimulation();
    prob.writeOutput("Couette.vtk");

    return 0;
}
```


