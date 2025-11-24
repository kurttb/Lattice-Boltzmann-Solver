#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "D2Q9Problem.hpp"
#include "VtkWriter.hpp"
#include "ComputeState.hpp"
#include "Collisions.hpp"
#include "Streaming.hpp"
#include "BoundaryConditions.hpp"

namespace LBM {

	/////////// Preprocessor Functions
	
	// Constructor
	D2Q9Problem::D2Q9Problem(const size_t Nx, const size_t Ny) {

		// Set Grid sizes
		_gridObj.Nx = Nx;
		_gridObj.Ny = Ny;

		// Set state sizes
		_rho.resize(Nx*Ny);
		_ux.resize(Nx*Ny);
		_uy.resize(Nx*Ny);

		// Set distribution function size
		_f.resize(9*Nx*Ny);

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
	void D2Q9Problem::setViscosity(const double nu) {
		_nu = nu;
	}


	// Set initial condition
	void D2Q9Problem::setIC(const double rho0, const double ux0, const double uy0) {

		// Fill vectors
		std::fill(_rho.begin(), _rho.end(), rho0);
		std::fill(_ux.begin(), _ux.end(), ux0);
		std::fill(_uy.begin(), _uy.end(), uy0);

		// Set distribution function to equilibrium
		for(size_t j = 0; j < _gridObj.Ny; ++j) {
			for(size_t i = 0; i < _gridObj.Nx; ++i) {
				size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
				double uxn = _ux[n]; // Extract x velocity component
				double uyn = _uy[n]; // Extract y velocity component
				double u_sq_ind = _ux[n]*_ux[n] + _uy[n]*_uy[n]; // Square velocity at the current grid node
				size_t base = 9*n;

				for(size_t k = 0; k < 9; ++k) {
					double e_dot_u = uxn*static_cast<double>(_ex[k]) + uyn*static_cast<double>(_ey[k]);
					_f[base + k] = _w[k] * _rho[n] * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
				}
			}
		}
	}

	// Set Boundary Conditions
	void D2Q9Problem::setBC(const string& BCName, const string& BCType, double uT) {
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
	void D2Q9Problem::setForces(const double Fx, const double Fy) {
		_Fx = Fx;
		_Fy = Fy;
	}






	// Solve
	void D2Q9Problem::runSimulation() {
				// Define knobs
		size_t N = _gridObj.Nx*_gridObj.Ny; // Total number of grid nodes


		// Derive characteristics of the flow physics
		constexpr double cs2 = 1.0 / 3.0; // Speed of sound squared
		double tau = ( _nu/(cs2) ) + 0.5; // Relaxation parameter
		double omega = 1.0 / tau;


		// Define Lattice - Start at rest, go east, and move counterclockwise
		//constexpr int ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
		//constexpr int ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
		//constexpr double w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distro



		// Allocate Distribution Functions
		vector<double> fstream(N*9, 0.0); // Distribution function after advection(9 for each grid node)


		// Initialize to equilibrium
	//	for(size_t j = 0; j < _gridObj.Ny; ++j) {
	//		for(size_t i = 0; i < _gridObj.Nx; ++i) {
	//			size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
	//			double uxn = _ux[n]; // Extract x velocity component
	//			double uyn = _uy[n]; // Extract y velocity component
	//			double u_sq_ind = _ux[n]*_ux[n] + _uy[n]*_uy[n]; // Square velocity at the current grid node
	//			size_t base = 9*n;

	//			for(size_t k = 0; k < 9; ++k) {
	//				double e_dot_u = uxn*static_cast<double>(_ex[k]) + uyn*static_cast<double>(_ey[k]);
	//				_f[base + k] = _w[k] * _rho[n] * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
	//			}
	//		}
	//	}


		// Boundary Indices
		//constexpr size_t ixL = 0;
		//size_t ixR = _gridObj.Nx - 1;
		//constexpr size_t iyB = 0;
		//size_t iyT = _gridObj.Ny - 1;


		// Start Update Loop

#pragma omp parallel default(shared)
	{
		for (size_t it = 0; it < _Nt; ++it) {

			// Compute macroscopic quantities from the distribution 
			D2Q9ReconstructState(_rho, _ux, _uy, _f, _ex, _ey, _gridObj, _Fx, _Fy, tau);


			// Collision step
			D2Q9BGKCollision(_rho, _ux, _uy, _f, _ex, _ey, _w, _gridObj, omega);


			// Streaming step 
			D2Q9Stream(_f, fstream, _ex, _ey, _gridObj);
			


			// Advection step on domain (Excluding top and right boundaries)
		//	#pragma omp for collapse(2) schedule(static)
		//	for(size_t j = 1; j < _gridObj.Ny - 1; ++j) {
		//		for(size_t i = 1; i < _gridObj.Nx - 1; ++i) {
		//			size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
		//			size_t base = 9*n;

		//			for(size_t k = 0; k < 9; ++k) {
		//				size_t i_new = i + ex[k];
		//				size_t j_new = j + ey[k];
		//				size_t n_new = i_new + _gridObj.Nx*j_new; // Location of point to be advected to

		//				fstream[9*n_new + k] = f[base + k];

		//			}
		//		}
		//	}


		//	// Streaming left/right
		//	#pragma omp for
		//	for(size_t j = 1; j < _gridObj.Ny - 1; ++j) {
		//		for(auto i: { size_t(0), size_t(_gridObj.Nx - 1) }) {
		//			size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
		//			size_t base = 9*n;

		//			for(size_t k = 0; k < 9; ++k) {
		//				size_t i_new = (i + ex[k] + _gridObj.Nx) % _gridObj.Nx;
		//				size_t j_new = j + ey[k];
		//				size_t n_new = i_new + _gridObj.Nx*j_new; // Location of point to be advected to

		//				fstream[9*n_new + k] = f[base + k];

		//			}
		//		}
		//	}


		//	// Streaming top/bottom
		//	#pragma omp for
		//	for(size_t i = 1; i < _gridObj.Nx - 1; ++i) {
		//		for(auto j: { size_t(0), size_t(_gridObj.Ny - 1) }) {
		//			size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
		//			size_t base = 9*n;

		//			for(size_t k = 0; k < 9; ++k) {
		//				size_t i_new = i + ex[k];
		//				size_t j_new = (j + ey[k] + _gridObj.Ny) % _gridObj.Ny;

		//				size_t n_new = i_new + _gridObj.Nx*j_new; // Location of point to be advected to

		//				fstream[9*n_new + k] = f[base + k];

		//			}
		//		}
		//	}



		//	// Stream the corners
		//	#pragma omp single
		//	{
		//		for(auto i: {size_t(0), size_t(_gridObj.Nx - 1)}) {
		//			for(auto j: { size_t(0), size_t(_gridObj.Ny - 1) }) {
		//				size_t n = i + _gridObj.Nx*j; // Current position on the flattened grid
		//				size_t base = 9*n;

		//				for(size_t k = 0; k < 9; ++k) {
		//					size_t i_new = (i + ex[k] + _gridObj.Nx) % _gridObj.Nx;
		//					size_t j_new = (j + ey[k] + _gridObj.Ny) % _gridObj.Ny;
		//					size_t n_new = i_new + _gridObj.Nx*j_new; // Location of point to be advected to

		//					fstream[9*n_new + k] = f[base + k];

		//				}
		//			}
		//		}

		//	}



			// Swap with f with fstream
			#pragma omp single
			{
				_f.swap(fstream);
			}


			// Enforce Boundary Conditions
			//enforceD2Q9BCs(_f, _w, _gridObj, U_lid, cs2, iyT); //Legacy

			if (_BCTop.BCType == "WallTangentVelocity") {
				tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCTop);
			}
			else if (_BCTop.BCType == "BounceBack") {
				bounceBackD2Q9(_f, _gridObj, _BCTop);
			}

			if (_BCBottom.BCType == "WallTangentVelocity") {
				tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCBottom);
			}
			else if (_BCBottom.BCType == "BounceBack") {
				bounceBackD2Q9(_f, _gridObj, _BCBottom);
			}

			if (_BCRight.BCType == "WallTangentVelocity") {
				tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCRight);
			}
			else if (_BCRight.BCType == "BounceBack") {
				bounceBackD2Q9(_f, _gridObj, _BCRight);
			}

			if (_BCLeft.BCType == "WallTangentVelocity") {
				tangentVelocityD2Q9(_f, _w, cs2, _gridObj, _BCLeft);
			}
			else if (_BCLeft.BCType == "BounceBack") {
				bounceBackD2Q9(_f, _gridObj, _BCLeft);
			}


			// Bounce-back on bottom wall
//			#pragma omp for
//			for (size_t i = 0; i < _gridObj.Nx; ++i) {
//
//				size_t n = i; // Current position on the flattened grid (j = 0)
//
//				f[9*n + 2] = f[9*n + 6];
//				f[9*n + 3] = f[9*n + 7];
//				f[9*n + 4] = f[9*n + 8];
//			}
//
//
//			// Moving lid on top wall
//			#pragma omp for
//			for (size_t i = 0; i < _gridObj.Nx; ++i) {
//				size_t n = i + _gridObj.Nx*iyT; // Current position on the flattened grid
//
//				// Calculate density at current grid location
//				double rho_ij = 0.0;
//				for (size_t k = 0; k < 9; ++k) {
//					double f_curr = f[9*n + k];
//					rho_ij += f_curr;
//
//				}
//
//				// Update distribution function
//				f[9*n + 6] = f[9*n + 2] - 2*w[6]*rho_ij*U_lid/cs2;
//				f[9*n + 7] = f[9*n + 3];
//				f[9*n + 8] = f[9*n + 4] + 2*w[8]*rho_ij*U_lid/cs2;
//			}
//
		}
	}

	}






	/////////// Post-processing functions
	// VTK Write
	void D2Q9Problem::writeOutput(string filePath) {
		std::string pv_title = "LBM Field";
		WriteVtk(_rho, _ux, _uy, _gridObj.Nx, _gridObj.Ny, filePath, pv_title);
	}



	// Destructor
	D2Q9Problem::~D2Q9Problem() = default;


}







