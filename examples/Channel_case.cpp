// Example case using the modularized library for the case of Channel flow
#include <string>
#include "D2Q9Problem.hpp"


int main() {

	// Grid Points in Each Direction
	const size_t Nx = 800;
	const size_t Ny = 100;

	// Set Problem Parameter
	const double Re = 100; // Reynolds Number
	const double Ma = 0.1; // Mach Number
	const size_t LChar = Ny; // Characteristic length scale

	// Derive Viscosity and Derive Characteristic Velocity
	double cs = 1.0 / sqrt(3); // Speed of sound
	double uChar = Ma*cs; // Characteristic Velocity
	double nu = uChar * LChar / Re; // Viscosity

	// Define D2Q9 Lattice Problem
	auto prob = LBM::D2Q9Problem(Nx, Ny);

	// Set Viscosity
	prob.setViscosity(nu); 

	// Initial Conditions
	const double rho0 = 1.0;
	const double ux0 = 0.0;
	const double uy0 = 0.0;
	prob.setIC(rho0, ux0, uy0);

	// Set time step
	size_t Nt = 50000;
	prob.setNumTimeSteps(Nt);

	// Set Boundary Conditions
	std::string BCTop = "Top"; // Options: Top, Bottom, Right, Left
	std::string BCTopType = "BounceBack"; // Options: Periodic, WallTangentVelocity, BounceBack
	prob.setBC(BCTop, BCTopType); // Sets the boundary condition

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
	const double Fx = 1e-5;
	const double Fy = 0.0;
	prob.setForces(Fx, Fy);


	// Run the simulation
	prob.runSimulation();

	// Write the output
	std::string filePath = "Channel.vtk";
	prob.writeOutput(filePath);


	return 0;
}
