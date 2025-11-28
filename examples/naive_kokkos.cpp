// naive_implementation.cpp with Kokkos
// 11/16/25
// A first transcription of the MATLAB Lattice Boltmzann Implemention for low-mach incompressible flows to C++

#include <Kokkos_Core.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "VtkWriter.hpp"

using namespace std;
using vec1 = Kokkos::View<double*>;
using vec2 =  Kokkos::View<double**>;
using vec3 =  Kokkos::View<double***>;


int main() {

	// Initialize Kokkos
	Kokkos::initialize();
	{
		// Define Knobs
		constexpr size_t Nx = 1000; // Number of x coordinates
		constexpr size_t Ny = 1000; // Number of y coordinates
		constexpr size_t N = Nx*Ny; // Total number of grid nodes
		constexpr double Ma = 0.1; // Mach number
		constexpr double Re = 100; // Reynolds number
		constexpr size_t max_it = 1000; // Maximum number of iterations

		// Initial Condition
		constexpr double rho_init = 1.0; // Initial density field
		constexpr double ux_init = 0.0; // Initial x velocity
		constexpr double uy_init = 0.0; // Initial y velocity

		// Derive characteristics of the flow physics
		constexpr size_t L = Ny - 1; // Length of the domain in the lattice
		const double cs = 1.0 / sqrt(3); // Speed of sound
		constexpr double cs2 = 1.0 / 3.0; // Speed of sound squared
		const double U_lid = Ma*cs; // Lid velocity
		const double nu = (U_lid*L) / Re; // Kinematic viscosity
		const double tau = ( nu/(cs2) ) + 0.5; // Relaxation parameter
		const double omega = 1.0 / tau;

		// Set body force
		double Fx = 0;

		// Define Lattice - Start at rest, go east, and move counterclockwise
		constexpr int ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
		constexpr int ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
		constexpr double w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distro
		//constexpr Kokkos::Array<int, 9> ex = {0, 1, 1, 0, -1, -1, -1, 0, 1};
		//constexpr Kokkos::Array<int, 9> ey = {0, 0, 1, 1, 1, 0, -1, -1, -1};
		//constexpr Kokkos::Array<double, 9> w = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distro

		// Allocate Macroscopic Fields
		Kokkos::View<double**> rho("rho", Ny, Nx);
		Kokkos::View<double**> ux("ux", Ny, Nx);
		Kokkos::View<double**> uy("uy", Ny, Nx);

		// Host device copies
		auto rho_host = Kokkos::create_mirror_view(rho);
		auto ux_host = Kokkos::create_mirror_view(ux);
		auto uy_host = Kokkos::create_mirror_view(uy);

		// Initalize Macroscopic Fields
		Kokkos::deep_copy(rho, rho_init);
		Kokkos::deep_copy(ux, ux_init);
		Kokkos::deep_copy(uy, uy_init);

		// Allocate Distribution Functions
		Kokkos::View<double***> f("f", Ny, Nx, 9); // Distribution funcs
		Kokkos::View<double***> fstream("fstream", Ny, Nx, 9); // Streaming distribution funcs


		// Initialize to equilibrium
		Kokkos::parallel_for("equilibrium_calc",
			Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {Ny, Nx}), 
			KOKKOS_LAMBDA(const int j, const int i) {
                double uxn = ux(j, i); // Extract x velocity component
                double uyn = uy(j, i); // Extract y velocity component
                double u_sq_ind = uxn*uxn + uyn*uyn; // Square velocity at the current grid

				for(int k = 0; k < 9; ++k) {
					double e_dot_u = uxn*static_cast<double>(ex[k]) + uyn*static_cast<double>(ey[k]);
					f(j, i, k) = w[k] * rho(j, i) * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);

				}
			}
        );

		
		// Set Boundary Indices
		constexpr size_t ixL = 0;
		constexpr size_t ixR = Nx - 1;
		constexpr size_t iyB = 0;
		constexpr size_t iyT = Ny - 1;

		
		// Loop over all time steps
		for (size_t it = 0; it < max_it; ++it) {

			// Compute macroscopic quantities from distribution
			Kokkos::parallel_for("compute_macroscopic", 
				Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {Ny, Nx}), 
				KOKKOS_LAMBDA(const int j, const int i) {
					double rho_ij = 0;
					double ux_ij = 0;
					double uy_ij = 0;

					for (int k = 0; k < 9; ++k) {
						double f_curr = f(j, i, k);
						rho_ij += f_curr;
						ux_ij += f_curr * ex[k];
						uy_ij += f_curr * ey[k];
					}


					// Update vals
					ux(j, i) = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
					uy(j, i) = uy_ij / rho_ij;
					rho(j, i) = rho_ij;
				}

			);


			// Collision Step
			Kokkos::parallel_for("compute_collision",
				Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {Ny, Nx}),
				KOKKOS_LAMBDA(const int j, const int i) {
					double uxn = ux(j, i);
					double uyn = uy(j, i);
					double u_sq_ind = uxn*uxn + uyn*uyn; // Square velocity at current grid node

					for (int k = 0; k < 9; ++k) {
						double e_dot_u = uxn*static_cast<double>(ex[k]) + uyn*static_cast<double>(ey[k]);
						double f_curr = f(j, i, k);
						double feq_curr = w[k] * rho(j, i) * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
						f(j, i, k) = f_curr - omega * (f_curr - feq_curr);
					}
				}
			);

			// Streaming Step on Interior
			Kokkos::parallel_for("Streaming_interior",
				Kokkos::MDRangePolicy<Kokkos::Rank<2>> ({1, 1},{Ny - 1, Nx - 1}),
				KOKKOS_LAMBDA(const int j, const int i) {
					for (int k = 0; k < 9; ++k) {
						int i_new = i + ex[k];
						int j_new = j + ey[k];
						fstream(j_new, i_new, k) = f(j, i, k);
					}
				}
			);

            // Streaming step left/right
			Kokkos::parallel_for("stream_lr",
				Kokkos::RangePolicy<>(1, Ny - 1),
				KOKKOS_LAMBDA(const int j) {

                    const int Nx_i = static_cast<int>(Nx);
					int i_vals[2] = {0, Nx_i - 1};

					for(int index = 0; index < 2; ++index) {
						int i = i_vals[index];
						for(int k = 0; k < 9; ++k) {
							int i_new = (i + ex[k] + Nx_i) % Nx_i;
							int j_new = j + ey[k];
							
							fstream(j_new, i_new, k) = f(j, i, k);
						}
					}
				}
			);


			// Streaming Top/Bottom
			Kokkos::parallel_for("stream_tb",
				Kokkos::RangePolicy<>(1, Nx - 1),
				KOKKOS_LAMBDA(const int i) {

					const int Ny_i = static_cast<int>(Ny);
					int j_vals[2] = {0, Ny_i - 1};

					for(int index = 0; index < 2; ++index) {
						int j = j_vals[index];
						for(int k = 0; k < 9; ++k) {
							int i_new = i + ex[k];
							int j_new = (j + ey[k] + Ny_i) % Ny_i;

							fstream(j_new, i_new, k) = f(j, i, k);
						}
					}
				}
			);

			// Stream the corners
			Kokkos::Array<int, 2> i_indices = {0, static_cast<int>(Nx - 1)};
			Kokkos::Array<int, 2> j_indices = {0, static_cast<int>(Ny - 1)};
			Kokkos::parallel_for("stream_corners",
				Kokkos::MDRangePolicy<Kokkos::Rank<2>> ({0, 0}, {2, 2}),
				KOKKOS_LAMBDA(const int i_index, const int j_index) {
					int i = i_indices[i_index];
				    int j = j_indices[j_index];
					int Nx_i = static_cast<int>(Nx);
					int Ny_i = static_cast<int>(Ny);

					for(size_t k = 0; k < 9; ++k) {
						int i_new = (i + ex[k] + Nx_i) % Nx_i;
						int j_new = (j + ey[k] + Ny_i) % Ny_i;

						fstream(j_new, i_new, k) = f(j, i, k);
					}
				}
			);

			// Swap fstream and f
            std::swap(f, fstream);

			// Bounce-back on bottom wall
			Kokkos::parallel_for("BounceBack_bot",
				static_cast<int>(Nx),
				KOKKOS_LAMBDA(const int i) {
					int j = 0;

					f(j, i, 2) = f(j, i, 6);
					f(j, i, 3) = f(j, i, 7);
					f(j, i, 4) = f(j, i, 8);
				}
			);


			// Moving lid on top wall
			Kokkos::parallel_for("Moving_top",
				static_cast<int>(Nx),
				KOKKOS_LAMBDA(const int i) {
					int j = Ny - 1;

					// Calculate density at current grid location
					double rho_ij = 0.0;
					for (int k = 0; k < 9; ++k) {
						double f_curr = f(j, i, k);
						rho_ij += f_curr;

					}

					// Update distribution function
					f(j, i, 6) = f(j, i, 2) - 2*w[6]*rho_ij*U_lid/cs2;
					f(j, i, 7) = f(j, i, 3);
					f(j, i, 8) = f(j, i, 4) + 2*w[8]*rho_ij*U_lid/cs2;
				}
			);
		// End of time loop
		}

		// Copy from device to host
		Kokkos::deep_copy(rho_host, rho);
		Kokkos::deep_copy(ux_host, ux);
		Kokkos::deep_copy(uy_host, uy);

		// Create vectors to load views into
		vector<double> rho_vec(Nx*Ny);
		vector<double> ux_vec(Nx*Ny);
		vector<double> uy_vec(Nx*Ny);

		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				rho_vec[i + Nx*j] = rho_host(j, i);
				ux_vec[i + Nx*j] = ux_host(j, i);
				uy_vec[i + Nx*j] = uy_host(j, i);
			}
		}


		// Write to file
		std::string file_path = "kokkos_couette.vtk";
		std::string pv_title = "LBM Field";
		LBM::WriteVtk(rho_vec, ux_vec, uy_vec, Nx, Ny, file_path, pv_title);

	}
	Kokkos::finalize();


	return 0;
}
