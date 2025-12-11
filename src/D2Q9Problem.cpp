#include <Kokkos_Core.hpp>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "GeneralProblem.hpp"
#include "D2Q9Problem.hpp"
#include "EquilibriumKokkos.hpp"
#include "ComputeStateKokkos.hpp"
#include "CollisionKokkos.hpp"
#include "StreamingKokkos.hpp"
#include "BCsKokkos.hpp"
#include "VtkWriter.hpp"

namespace LBM {

	/////////// Preprocessor Functions
	
	// Constructor
	D2Q9Problem::D2Q9Problem(const int Nx, const int Ny) : GeneralProblem()
	{

		if (Nx == 0 || Ny == 0) {
			throw std::invalid_argument("Grid dimensions must be non-zero");
		}

		// Set Grid sizes
		_gridObj.Nx = Nx;
		_gridObj.Ny = Ny;

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


	// Set Boundary Conditions
	void D2Q9Problem::setBC(const std::string& BCName, const std::string& BCType, float uT) {
		if (BCName == "Top") {
			_BCTop.BCType = BCType;
			_BCTop.U_wall = uT;
		}
		else if (BCName == "Bottom") {
			_BCBottom.BCType = BCType;
			_BCBottom.U_wall = uT;
		}
		else if (BCName == "Right") {
			_BCRight.BCType = BCType;
			_BCRight.U_wall = uT;
		}
		else if (BCName == "Left") {
			_BCLeft.BCType = BCType;
			_BCLeft.U_wall = uT;
		}
		else {
			cout << "Boundary Name is not recognized" << endl;
			exit(EXIT_FAILURE);
		}
	}


	// Solve
	void D2Q9Problem::runSimulation() {
		
		bool initialized_here = false;
		if (!Kokkos::is_initialized()) {
			Kokkos::initialize();
			initialized_here = true;
		}
		{
			// Show what Kokkos is running on
			Kokkos::print_configuration(std::cout);
			
			// Define knobs
			int Nx = _gridObj.Nx;
			int Ny = _gridObj.Ny;
			int N = Nx*Ny; // Total number of grid nodes
			
			// Derive characteristics of the flow physics
			const float cs2 = 1.0f / 3.0f; // Speed of sound squared
			const float tau = ( _nu/(cs2) ) + 0.5f; // Relaxation parameter
			const float omega = 1.0f / tau;


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
				CalcEq(rho, ux, uy, f)
			);


			// Set timer
			Kokkos::Timer timer;


			// Start Update Loop
			for (int it = 0; it < _Nt; ++it) {

				// Compute macroscopic quantities from the distribution 
				Kokkos::parallel_for("ComputeState",
					N,
					ComputeState(rho, ux, uy, f, _Fx, _Fy, tau)
				);
				

				// Collision step
				Kokkos::parallel_for("Collision",
					N,
					ComputeCollision(rho, ux, uy, f, omega)
				);


				// Streaming step 
				Kokkos::parallel_for("Streaming",
					N,
					ComputeStreaming(f, fstream, Nx, Ny)
				);
				


				// Swap with f with fstream
				std::swap(f, fstream);


				// Enforce Boundary Conditions
				if (_BCTop.BCType == "WallTangentVelocity") {
					Kokkos::parallel_for("TopTanVel",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCTop.j_min, _BCTop.i_min}, {_BCTop.j_max + 1, _BCTop.i_max + 1}),
						TangentVelocity(f, _BCTop.f_inc, _BCTop.f_ref, _BCTop.U_wall, cs2, Nx)
					);
				}
				else if (_BCTop.BCType == "BounceBack") {
					Kokkos::parallel_for("TopBounceBack",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCTop.j_min, _BCTop.i_min}, {_BCTop.j_max + 1, _BCTop.i_max + 1}),
						BounceBack(f, _BCTop.f_inc, _BCTop.f_ref, Nx)
					);
				}

				if (_BCBottom.BCType == "WallTangentVelocity") {
					Kokkos::parallel_for("BotTanVel",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCBottom.j_min, _BCBottom.i_min}, {_BCBottom.j_max + 1, _BCBottom.i_max + 1}),
						TangentVelocity(f, _BCBottom.f_inc, _BCBottom.f_ref, _BCBottom.U_wall, cs2, Nx)
					);
				}
				else if (_BCBottom.BCType == "BounceBack") {
					Kokkos::parallel_for("BottomBounceBack",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCBottom.j_min, _BCBottom.i_min}, {_BCBottom.j_max + 1, _BCBottom.i_max + 1}),
						BounceBack(f, _BCBottom.f_inc, _BCBottom.f_ref, Nx)
					);
				}

				if (_BCRight.BCType == "WallTangentVelocity") {
					Kokkos::parallel_for("RightTanVel",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCRight.j_min, _BCRight.i_min}, {_BCRight.j_max + 1, _BCRight.i_max + 1}),
						TangentVelocity(f, _BCRight.f_inc, _BCRight.f_ref, _BCRight.U_wall, cs2, Nx)
					);
				}
				else if (_BCRight.BCType == "BounceBack") {
					Kokkos::parallel_for("RightBounceBack",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCRight.j_min, _BCRight.i_min}, {_BCRight.j_max + 1, _BCRight.i_max + 1}),
						BounceBack(f, _BCRight.f_inc, _BCRight.f_ref, Nx)
					);
				}

				if (_BCLeft.BCType == "WallTangentVelocity") {
					Kokkos::parallel_for("LeftTanVel",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCLeft.j_min, _BCLeft.i_min}, {_BCLeft.j_max + 1, _BCLeft.i_max + 1}),
						TangentVelocity(f, _BCLeft.f_inc, _BCLeft.f_ref, _BCLeft.U_wall, cs2, Nx)
					);
				}
				else if (_BCLeft.BCType == "BounceBack") {
					Kokkos::parallel_for("LeftBounceBack",
						Kokkos::MDRangePolicy<Kokkos::Rank<2>>({_BCLeft.j_min, _BCLeft.i_min}, {_BCLeft.j_max + 1, _BCLeft.i_max + 1}),
						BounceBack(f, _BCLeft.f_inc, _BCLeft.f_ref, Nx)
					);
				}

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
		if (initialized_here) {
			Kokkos::finalize();
		}
	}






	/////////// Post-processing functions
	// VTK Write
	void D2Q9Problem::writeOutput(std::string filePath) {
		std::string pv_title = "LBM Field";
		WriteVtk(_rho, _ux, _uy, _gridObj.Nx, _gridObj.Ny, filePath, pv_title);
	}

	// Destructor
	D2Q9Problem::~D2Q9Problem() = default;

	// Getters
	std::vector<float> D2Q9Problem::getRho() const { return _rho; }
	std::vector<float> D2Q9Problem::getUx() const { return _ux; }
	std::vector<float> D2Q9Problem::getUy() const { return _uy; }

}







