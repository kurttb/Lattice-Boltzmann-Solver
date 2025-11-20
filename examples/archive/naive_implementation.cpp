// naive_implementation.cpp
// 11/9/25
// A first transcription of the MATLAB Lattice Boltmzann Implemention for low-mach incompressible flows to C++

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "vtk_writer.hpp"

using namespace std;


int main() {

	// Define knobs
	string flow_type = "Couette"; 
	const size_t Nx = 100; // Number of x coordinates
	const size_t Ny = 100; // Number of y coordinates
	const size_t N = Nx*Ny; // Total number of grid nodes
	const double Ma = 0.1; // Mach number
	const double Re = 100; // Reynolds number
	const size_t max_it = 50000; // Maximum number of iterations
	const double tol = 1e-4; // Steady-state tolerance

	// Initial Condition
	double rho_init = 1.0; // Initial density field
	double u_init[2] = {0.0, 0.0}; // Initial Velocity

	// Derive characteristics of the flow physics
	const size_t L = Ny - 1; // Length of the domain in the lattice
	const double cs = 1.0 / sqrt(3); // Speed of sound
	const double cs2 = 1.0 / 3.0; // Speed of sound squared
	const double U_lid = Ma*cs; // Lid velocity
	const double nu = (U_lid*L) / Re; // Kinematic viscosity
	const double tau = ( nu/(cs2) ) + 0.5; // Relaxation parameter
	const double omega = 1.0 / tau;


	double Fx = 0;
	if (flow_type == "Channel") {
		Fx += 1e-4;
	}

	// Define Lattice - Start at rest, go east, and move counterclockwise
	const int ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
	const int ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
	const double w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distro


	// Allocate Macroscopic Fields
	vector<double> rho(N, rho_init); // Density field
	vector<double> ux(N, u_init[0]); // x-Velocity field
	vector<double> uy(N, u_init[1]); // y-Velocity field

	// Allocate Distribution Functions
	vector<double> f(N*9, 0.0); // Distribution functions (9 for each grid node)
	vector<double> fcoll(N*9, 0.0); // Distribution function after collision(9 for each grid node)
	vector<double> fstream(N*9, 0.0); // Distribution function after advection(9 for each grid node)


	// Calculate squared speed at each grid node
	vector<double> u_sq(N);

	for(size_t i = 0; i < N; ++i) {
		u_sq.at(i) = ux.at(i)*ux.at(i) + uy.at(i)*uy.at(i);
	}


	// Initialize to equilibrium
	for(size_t j = 0; j < Ny; ++j) {
		for(size_t i = 0; i < Nx; ++i) {
			size_t n = i + Nx*j; // Current position on the flattened grid
			double uxn = ux.at(n); // Extract x velocity component
			double uyn = uy.at(n); // Extract y velocity component
			double u_sq_ind = u_sq.at(n); // Square velocity at the current grid node
			for(size_t k = 0; k < 9; ++k) {
				double e_dot_u = uxn*static_cast<double>(ex[k]) + uyn*static_cast<double>(ey[k]);
				f.at(9*n + k) = w[k] * rho.at(n) * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
			}
		}
	}

/*	for (size_t i = 0; i < f.size(); i++ ) {
		cout << f.at(i) << endl;
	}
	exit(0); */

	// Boundary Indices
	size_t ixL = 0;
	size_t ixR = Nx - 1;
	size_t iyB = 0;
	size_t iyT = Ny - 1;


	// Start Update Loop
	vector<double> prev_ux = ux; // ux at the previous iteration
	size_t it = 0; // Current iteration

	while (it < max_it) {

		// Update iteration
		it += 1;

		// Compute macroscopic quantities from distribution
		for(size_t j = 0; j < Ny; ++j) {
			for (size_t i = 0; i < Nx; ++i) {
				size_t n = i + Nx*j; // Current position on the flattened grid
				double rho_ij = 0;
				double ux_ij = 0;
				double uy_ij = 0;

				for (size_t k = 0; k < 9; ++k) {
					double f_curr = f.at(9*n + k);
					rho_ij += f_curr;
					ux_ij += f_curr * ex[k];
					uy_ij += f_curr * ey[k];
				}

				ux.at(n) = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
				uy.at(n) = uy_ij / rho_ij;
				rho.at(n) = rho_ij;
			}
		}

		/* for (size_t i = 0; i < ux.size(); ++i) {
			cout << rho.at(i) << endl;
		}

		exit(0); */


		// Recompute square velocity
		for(size_t i = 0; i < N; ++i) {
			u_sq.at(i) = ux.at(i)*ux.at(i) + uy.at(i)*uy.at(i);
		}

		// Collision Step
		for(size_t j = 0; j < Ny; ++j) {
			for(size_t i = 0; i < Nx; ++i) {
				size_t n = i + Nx*j; // Current position on the flattened grid
				double uxn = ux.at(n); // Extract x velocity component
				double uyn = uy.at(n); // Extract y velocity component
				double u_sq_ind = u_sq.at(n); // Square velocity at the current grid node
				for(size_t k = 0; k < 9; ++k) {
					double e_dot_u = uxn*static_cast<double>(ex[k]) + uyn*static_cast<double>(ey[k]);
					double f_curr = f.at(9*n + k);
					double feq_curr = w[k] * rho.at(n) * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind); // Calculate equilbrium distribution
					f.at(9*n + k) = f_curr - omega * (f_curr - feq_curr);
				}
			}
		}



		/*for (size_t i = 0; i < f.size(); i++ ) {
			cout << f.at(i) << endl;
		}
		exit(0); */


		// Advection step on domain (Excluding top and right boundaries)
		fill(fstream.begin(), fstream.end(), 0.0);
		for(size_t j = 1; j < Ny - 1; ++j) {
			for(size_t i = 1; i < Nx - 1; ++i) {
				size_t n = i + Nx*j; // Current position on the flattened grid

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = i + ex[k];
					size_t j_new = j + ey[k];
					size_t n_new = i_new + Nx*j_new; // Location of point to be advected to

					fstream.at(9*n_new + k) = f.at(9*n + k);

				}
			}
		}


		/* for (size_t i = 0; i < fstream.size(); i++ ) {
			cout << fstream.at(i) << endl;
		}
		exit(0); */


		// Streaming left/right
		for(size_t j = 1; j < Ny - 1; ++j) {
			for(auto i: { size_t(0), size_t(Nx - 1) }) {
				size_t n = i + Nx*j; // Current position on the flattened grid

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = (i + ex[k] + Nx) % Nx;
					size_t j_new = j + ey[k];
					size_t n_new = i_new + Nx*j_new; // Location of point to be advected to

					fstream.at(9*n_new + k) = f.at(9*n + k);

				}
			}
		}


		// Streaming top/bottom
		for(size_t i = 1; i < Nx - 1; ++i) {
			for(auto j: { size_t(0), size_t(Ny - 1) }) {
				size_t n = i + Nx*j; // Current position on the flattened grid

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = i + ex[k];
					size_t j_new = (j + ey[k] + Ny) % Ny;

					size_t n_new = i_new + Nx*j_new; // Location of point to be advected to

					fstream.at(9*n_new + k) = f.at(9*n + k);

				}
			}
		}

		 /*for (size_t i = 0; i < fstream.size(); i++ ) {
			cout << fstream.at(i) << endl;
		}
		exit(0); */


		// Stream the corners
		
		for(auto i: {size_t(0), size_t(Nx - 1)}) {
			for(auto j: { size_t(0), size_t(Ny - 1) }) {
				size_t n = i + Nx*j; // Current position on the flattened grid

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = (i + ex[k] + Nx) % Nx;
					size_t j_new = (j + ey[k] + Ny) % Ny;
					size_t n_new = i_new + Nx*j_new; // Location of point to be advected to

					fstream.at(9*n_new + k) = f.at(9*n + k);

				}
			}
		} 


		 /*for (size_t i = 0; i < fstream.size(); i++ ) {
			cout << fstream.at(i) << endl;
		 }
		 exit(0); */

		// Swap with f with fstream
		f.swap(fstream);

		// Bounce-back on bottom wall
		for (size_t i = 0; i < Nx; ++i) {

			size_t n = i; // Current position on the flattened grid (j = 0)

			f.at(9*n + 2) = f.at(9*n + 6);
			f.at(9*n + 3) = f.at(9*n + 7);
			f.at(9*n + 4) = f.at(9*n + 8);
		}


		// Moving lid on top wall
		for (size_t i = 0; i < Nx; ++i) {
			size_t n = i + Nx*iyT; // Current position on the flattened grid

			// Calculate density at current grid location
			double rho_ij = 0.0;
			for (size_t k = 0; k < 9; ++k) {
				double f_curr = f.at(9*n + k);
				rho_ij += f_curr;

			}

			// Update distribution function
			f.at(9*n + 6) = f.at(9*n + 2) - 2*w[6]*rho_ij*U_lid/cs2;
			f.at(9*n + 7) = f.at(9*n + 3);
			f.at(9*n + 8) = f.at(9*n + 4) + 2*w[8]*rho_ij*U_lid/cs2;
		}

	}

	// Print out x velocity values
	//for (size_t i = 0; i < ux.size(); i++) {
	//	cout << ux[i] << endl;
	//}

	// Write output
	std::string file_path = "tmp/paraview_ser.vtk";
	std::string pv_title = "LBM Field";
	write_vtk(rho, ux, uy, Nx, Ny, file_path, pv_title);






	return 0;

	


}
