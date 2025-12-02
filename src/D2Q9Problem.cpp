#include <Kokkos_Core.hpp>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "D2Q9Problem.hpp"
#include "EquilibriumKokkos.hpp"
#include "ComputeStateKokkos.hpp"
#include "CollisionKokkos.hpp"
#include "StreamingKokkos.hpp"
#include "VtkWriter.hpp"
#include "BoundaryConditions.hpp"

namespace LBM {

	/////////// Preprocessor Functions
	
	// Constructor
	D2Q9Problem::D2Q9Problem(const size_t Nx, const size_t Ny) {

		// Set Grid sizes
		_gridObj.Nx = Nx;
		_gridObj.Ny = Ny;

		// Set state sizes
		//_rho.resize(Nx*Ny);
		//_ux.resize(Nx*Ny);
		//_uy.resize(Nx*Ny);

		// Set distribution function size
		//_f.resize(9*Nx*Ny);

		// Set up top boundary condition
		_BCTop.f_inc[0] = 2;
		_BCTop.f_inc[1] = 3;
		_BCTop.f_inc[2] = 4;
		_BCTop.f_ref[0] = 6;
		_BCTop.f_ref[1] = 7;
		_BCTop.f_ref[2] = 8;

		_BCTop.i_min = 0;
		_BCTop.i_max = Nx - 1;
		_BCTop.j_min = Ny - 1;
		_BCTop.j_max = Ny - 1;


		// Set up Bottom boundary condition
		_BCBottom.f_inc[0] = 8;
		_BCBottom.f_inc[1] = 7;
		_BCBottom.f_inc[2] = 6;
		_BCBottom.f_ref[0] = 4;
		_BCBottom.f_ref[1] = 3;
		_BCBottom.f_ref[2] = 2;

		_BCBottom.i_min = 0;
		_BCBottom.i_max = Nx - 1;
		_BCBottom.j_min = 0;
		_BCBottom.j_max = 0;

		// Set up Right boundary condition
		_BCRight.f_inc[0] = 2;
		_BCRight.f_inc[1] = 1;
		_BCRight.f_inc[2] = 8;
		_BCRight.f_ref[0] = 6;
		_BCRight.f_ref[1] = 5;
		_BCRight.f_ref[2] = 4;


		_BCRight.i_min = Nx - 1;
		_BCRight.i_max = Nx - 1;
		_BCRight.j_min = 0;
		_BCRight.j_max = Ny - 1;

		// Set up Left boundary condition
		_BCLeft.f_inc[0] = 4;
		_BCLeft.f_inc[1] = 5;
		_BCLeft.f_inc[2] = 6;
		_BCLeft.f_ref[0] = 8;
		_BCLeft.f_ref[1] = 1;
		_BCLeft.f_ref[2] = 2;

		_BCLeft.i_min = 0;
		_BCLeft.i_max = 0;
		_BCLeft.j_min = 0;
		_BCLeft.j_max = Ny - 1;
	}

	// Set Viscosity
	void D2Q9Problem::setViscosity(const float nu) {
		_nu = nu;
	}


	// Set initial condition
	void D2Q9Problem::setIC(const float rho0, const float ux0, const float uy0) {
		_rho0 = rho0;
		_ux0 = ux0;
		_uy0 = uy0;

		// Fill vectors
		//std::fill(_rho.begin(), _rho.end(), rho0);
		//std::fill(_ux.begin(), _ux.end(), ux0);
		//std::fill(_uy.begin(), _uy.end(), uy0);

		// Set distribution function to equilibrium
		//for(size_t j = 0; j < _gridObj.Ny; ++j) {
			//for(size_t i = 0; i < _gridObj.Nx; ++i) {
				//size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
				//double uxn = _ux[n]; // Extract x velocity component
				//double uyn = _uy[n]; // Extract y velocity component
				//double u_sq_ind = _ux[n]*_ux[n] + _uy[n]*_uy[n]; // Square velocity at the current grid node
				//size_t base = 9*n;

				//for(size_t k = 0; k < 9; ++k) {
				//	double e_dot_u = uxn*static_cast<double>(_ex[k]) + uyn*static_cast<double>(_ey[k]);
				//	_f[base + k] = _w[k] * _rho[n] * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
				//}
			//}
		//}
	}

	// Set Boundary Conditions
	void D2Q9Problem::setBC(const std::string& BCName, const std::string& BCType, float uT) {
		if (BCName == "Top") {
			_BCTop.BCType = BCType;
			_BCTop.U_wall = uT;
			//cout << "Hello from the Top" << endl;
			//cout << _BCTop.U_wall << endl;
		}
		else if (BCName == "Bottom") {
			_BCBottom.BCType = BCType;
			_BCBottom.U_wall = uT;
			//cout << "Hello from the Bottom" << endl;
			//cout << _BCBottom.U_wall << endl;
		}
		else if (BCName == "Right") {
			_BCRight.BCType = BCType;
			_BCRight.U_wall = uT;
			//cout << "Hello from the Right" << endl;
			//cout << _BCRight.U_wall << endl;
		}
		else if (BCName == "Left") {
			_BCLeft.BCType = BCType;
			_BCLeft.U_wall = uT;
			//cout << "Hello from the Left" << endl;
			//cout << _BCLeft.U_wall << endl;
		}
		else {
			cout << "Boundary Name is not recognized" << endl;
			exit(EXIT_FAILURE);
		}
	}




	// Set number of time steps
	void D2Q9Problem::setNumTimeSteps(const size_t Nt) {
		_Nt = Nt;
	}

	// Set forces
	void D2Q9Problem::setForces(const float Fx, const float Fy) {
		_Fx = Fx;
		_Fy = Fy;
	}






	// Solve
	void D2Q9Problem::runSimulation() {
		
		Kokkos::initialize();
		{
			// Define knobs
			size_t Nx = _gridObj.Nx;
			size_t Ny = _gridObj.Ny;
			size_t N = Nx*Ny; // Total number of grid nodes

			// Forces
			const float Fx = _Fx;
			const float Fy = _Fy;

			// Derive characteristics of the flow physics
			const float U_lid = 0.058;
			const float cs2 = 1.0f / 3.0f; // Speed of sound squared
			const float tau = ( _nu/(cs2) ) + 0.5f; // Relaxation parameter
			const float omega = 1.0f / tau;

			// Allocate discrete velocity directions and weights
		//	Kokkos::View<int*> ex("ex", 9);
		//	Kokkos::View<int*> ey("ey", 9);
		//	Kokkos::View<float*> w("w", 9);
		//	auto ex_h = Kokkos::create_mirror_view(ex);
		//	auto ey_h = Kokkos::create_mirror_view(ey);
		//	auto w_h = Kokkos::create_mirror_view(w);

		//	// Fill host and device velocities/weights
		//	for (size_t i = 0; i < 9; ++i) {
		//		ex_h(i) = _ex[i];
		//		ey_h(i) = _ey[i];
		//	}

		//	Kokkos::deep_copy(ex, ex_h);
		//	Kokkos::deep_copy(ey, ey_h);
		//	Kokkos::deep_copy(w, w_h);

			constexpr int ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
			constexpr int ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
			constexpr float w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distribution
		//static const int* ex = _ex;
		//static const int* ey = _ey;
		//static const float* w = _w;


			// Allocate state and distribution function
			Kokkos::View<float*> rho("rho", N);
			Kokkos::View<float*> ux("ux", N);
			Kokkos::View<float*> uy("uy", N);
			Kokkos::View<float**> f("f", N, 9);

			// Initalize state and distribution function
			Kokkos::deep_copy(rho, _rho0);
			Kokkos::deep_copy(ux, _ux0);
			Kokkos::deep_copy(uy, _uy0);
			Kokkos::deep_copy(f, 0.0f);


			// Allocate and Initialize Streaming Distribution Function
			Kokkos::View<float**> fstream("fstream", N, 9);
			Kokkos::deep_copy(fstream, 0.0f);

			// Initialize Distribution Function to Equilibrium
			Kokkos::parallel_for("InitEq",
				N,
				KOKKOS_LAMBDA(const int n) {
					float uxn = ux(n);
					float uyn = uy(n);
					float rho_n = rho(n);
					float u_sq_ind = uxn*uxn + uyn*uyn;

					for (int k = 0; k < 9; ++k) {
						float e_dot_u = uxn*static_cast<float>(ex[k]) + uyn*static_cast<float>(ey[k]);
						f(n, k) = w[k] * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
					}
				}
			);
			
	//		Kokkos::parallel_for("InitEq",
	//			N,
	//			CalcEq(rho, ux, uy, f, _ex, _ey, _w)
	//		);

			// Boundary tag
			size_t iyT = Ny - 1;

			// Set timer
			Kokkos::Timer timer;


			// Start Update Loop
			for (size_t it = 0; it < _Nt; ++it) {

				// Compute macroscopic quantities from the distribution 
				//D2Q9ReconstructState(_rho, _ux, _uy, _f, _ex, _ey, _gridObj, _Fx, _Fy, tau);
				Kokkos::parallel_for("ComputeMacro",
					N,
					KOKKOS_LAMBDA(const int n) {
						float rho_ij = 0;
						float ux_ij = 0;
						float uy_ij = 0;

						for (int k = 0; k < 9; ++k) {
							float f_curr = f(n, k);
							rho_ij += f_curr;
							ux_ij += f_curr * ex[k];
							uy_ij += f_curr * ey[k];
						}

						ux(n) = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
						uy(n) = (uy_ij / rho_ij) + (Fy*tau / rho_ij);
						uy(n) = uy_ij / rho_ij;
						rho(n) = rho_ij;
					}
				);

		//		Kokkos::parallel_for("ComputeState",
		//			N,
		//			ComputeState(rho, ux, uy, f, _ex, _ey, _Fx, _Fy, tau)
		//		);
				


				// Collision step
				//D2Q9BGKCollision(_rho, _ux, _uy, _f, _ex, _ey, _w, _gridObj, omega);
			//	Kokkos::parallel_for("Collision",
			//		N,
			//		KOKKOS_LAMBDA(const int n) {
			//			float uxn = ux(n);
			//			float uyn = uy(n);
			//			float rho_n = rho(n);
			//			float u_sq_ind = uxn*uxn + uyn*uyn;

			//			for (int k = 0; k < 9; ++k) {
			//				float e_dot_u = uxn*static_cast<float>(ex[k]) + uyn*static_cast<float>(ey[k]);
			//				float f_curr = f(n, k);
			//				float feq_curr = w[k] * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
			//				f(n, k) = f_curr - omega * (f_curr - feq_curr);
			//			}
			//		}
			//	);

				Kokkos::parallel_for("Collision",
					N,
					ComputeCollision(rho, ux, uy, f, omega)
				);


				// Streaming step 
				//D2Q9Stream(_f, fstream, _ex, _ey, _gridObj);
				Kokkos::parallel_for("Streaming",
					N,
					KOKKOS_LAMBDA(const int n) {
						int i = n % Nx;
						int j = n / Nx;

						for (int k = 0; k < 9; ++k) {
							int i_dest = (i + ex[k] + Nx) % Nx;
							int j_dest = (j + ey[k] + Ny) % Ny;

							int n_new = i_dest + Nx*j_dest;

							fstream(n_new, k) = f(n, k);
						}
					}
				);

			//Kokkos::parallel_for("Streaming",
			//	N,
			//	ComputeStreaming(f, fstream, _ex, _ey, Nx, Ny)
			//);
				


				// Swap with f with fstream
				//#pragma omp single
				//{
				//	_f.swap(fstream);
				//}
				std::swap(f, fstream);


				// Enforce Boundary Conditions
			//	if (_BCTop.BCType == "WallTangentVelocity") {
			//		tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCTop);
			//	}
			//	else if (_BCTop.BCType == "BounceBack") {
			//		bounceBackD2Q9(_f, _gridObj, _BCTop);
		//			Kokkos::parallel_for("TopBounceBack",
		//				Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCTop.j_min, _BCTop.i_min}, {_BCTop.j_max + 1, _BCTop.i_max + 1}),
		//				BounceBack(f, _BCTop.f_inc, _BCTop.f_ref, Nx)
		//			);
			//	}

			//	if (_BCBottom.BCType == "WallTangentVelocity") {
			//		tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCBottom);
			//	}
			//	else if (_BCBottom.BCType == "BounceBack") {
			//		bounceBackD2Q9(_f, _gridObj, _BCBottom);
			//	}

			//	if (_BCRight.BCType == "WallTangentVelocity") {
			//		tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCRight);
			//	}
			//	else if (_BCRight.BCType == "BounceBack") {
			//		bounceBackD2Q9(_f, _gridObj, _BCRight);
			//	}

			//	if (_BCLeft.BCType == "WallTangentVelocity") {
			//		tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCLeft);
			//	}
			//	else if (_BCLeft.BCType == "BounceBack") {
			//		bounceBackD2Q9(_f, _gridObj, _BCLeft);
			//	}

				// BounceBack on bottom
				Kokkos::parallel_for("BounceBackBottom",
					Nx,
					KOKKOS_LAMBDA(const int i) {
						int n = i; // j = 0
						f(n, 2) = f(n, 6);
						f(n, 3) = f(n, 7);
						f(n, 4) = f(n, 8);
					}
				);

				// Moving lid on top
				Kokkos::parallel_for("MovingLidTop",
					Nx,
					KOKKOS_LAMBDA(const int i) {
						int n = i + Nx*iyT;
						float rho_ij = 0.0f;

						for (int k = 0; k < 9; ++k) {
							rho_ij += f(n, k);
						}
						f(n, 6) = f(n, 2) - 2.0f*w[6]*rho_ij*U_lid/cs2;
						f(n, 7) = f(n, 3);
						f(n, 8) = f(n, 4) + 2.0f*w[8]*rho_ij*U_lid/cs2;
					}
				);

			// End of time loop
			}

			Kokkos::fence();
			double time = timer.seconds();
			std::cout << "Time" << time << std::endl;

			// Fill vectors with the Kokkos views
			auto rho_h = Kokkos::create_mirror_view(rho);
			auto ux_h = Kokkos::create_mirror_view(ux);
			auto uy_h = Kokkos::create_mirror_view(uy);

			Kokkos::deep_copy(rho_h, rho);
			Kokkos::deep_copy(ux_h, ux);
			Kokkos::deep_copy(uy_h, uy);

			_rho.resize(Nx*Ny);
			_ux.resize(Nx*Ny);
			_uy.resize(Nx*Ny);
			for (int i = 0; i < N; ++i) {
				_rho[i] = rho_h(i);
				_ux[i] = ux_h(i);
				_uy[i] = uy_h(i);
			}

		}
		Kokkos::finalize();
	}






	/////////// Post-processing functions
	// VTK Write
	void D2Q9Problem::writeOutput(std::string filePath) {
		std::string pv_title = "LBM Field";
		WriteVtk(_rho, _ux, _uy, _gridObj.Nx, _gridObj.Ny, filePath, pv_title);
	}

	// Destructor
	D2Q9Problem::~D2Q9Problem() = default;


}







