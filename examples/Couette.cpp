// Example case using the modularized library for the case of Couette flow
#include <string>
#include "D2Q9Problem.hpp"


int main() {

	// # grid points in each direction
	const size_t Nx = 100;
	const size_t Ny = 100;

	// Define D2Q9 Lattice Problem
	auto prob = LBM::D2Q9Problem(Nx, Ny);

	// Initial conditons
	const double rho0 = 1.0;
	const double ux0 = 0.0;
	const double uy0 = 0.0;
	prob.setIC(rho0, ux0, uy0);

	// Set time step
	size_t Nt = 50000;
	prob.setNumTimeSteps(Nt);

	// Set body forces
	const double Fx = 0.0;
	const double Fy = 0.0;
	prob.setForces(Fx, Fy);


	// Run the simulation
	prob.runSimulation();

	// Write the output
	std::string filePath = "Couette.vtk";
	prob.writeOutput(filePath);


	return 0;
}
