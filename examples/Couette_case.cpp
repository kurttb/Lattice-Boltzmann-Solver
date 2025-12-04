// Example case using the modularized library for the case of Couette flow
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "D2Q9Problem.hpp"


int main(int argc, char* argv[]) {

	// Grid Points in Each Direction
	int Nx = 1000;
	int Ny = 1000;
	int Nt = 1000;

	if (argc >= 3) {
		Nx = std::atoi(argv[1]);
		Ny = std::atoi(argv[2]);
	}
	if (argc >= 4) {
		Nt = std::atoi(argv[3]);
	}

	std::cout << "Running with Nx: " << Nx << ", Ny: " << Ny << ", Nt: " << Nt << std::endl;

	// Set Problem Parameter
	const float Re = 100; // Reynolds Number
	const float Ma = 0.1; // Mach Number
	const int LChar = Ny; // Characteristic length scale

	// Derive Viscosity and Derive Characteristic Velocity
	float cs = 1.0 / std::sqrt(3); // Speed of sound
	float uChar = Ma*cs; // Characteristic Velocity
	float nu = uChar * LChar / Re; // Viscosity

	// Define D2Q9 Lattice Problem
	auto prob = LBM::D2Q9Problem(Nx, Ny);

	// Set Viscosity
	prob.setViscosity(nu); 

	// Initial Conditions
	const float rho0 = 1.0;
	const float ux0 = 0.0;
	const float uy0 = 0.0;
	prob.setIC(rho0, ux0, uy0);

	// Set time step
	prob.setNumTimeSteps(Nt);

	// Set Boundary Conditions
	std::string BCTop = "Top"; // Options: Top, Bottom, Right, Left
	std::string BCTopType = "WallTangentVelocity"; // Options: Periodic, WallTangentVelocity, BounceBack
	float uT = uChar; // Tangent Velocity
	prob.setBC(BCTop, BCTopType, uT); // Sets the boundary condition

	std::string BCBottom = "Bottom"; // Bottom BC
	std::string BCBottomType = "BounceBack"; // Bounce-back (equivalent to no-slip)
	prob.setBC(BCBottom, BCBottomType); // Set BC

	// NOTE: Periodic boundaries are implied if another BC is not specified. Periodic can be optionally specified
	std::string BCLeft = "Left";
	std::string BCLeftType = "Periodic";
	prob.setBC(BCLeft, BCLeftType);

	std::string BCRight = "Right";
	std::string BCRightType = "Periodic";
	prob.setBC(BCRight, BCRightType);
	

	// Set body forces
	const float Fx = 0.0;
	const float Fy = 0.0;
	prob.setForces(Fx, Fy);


	// Run the simulation
	prob.runSimulation();

	// Write the output
	std::string filePath = "Couette.vtk";
	prob.writeOutput(filePath);


	return 0;
}
