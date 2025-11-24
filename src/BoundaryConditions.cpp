// Implementation of boundary condition functions
#include <vector>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"

using namespace std;

namespace LBM {

	// Bounce-Back BC
	void bounceBackD2Q9(vector<double>& f,
						const LBM::CartesianGrid2D& gridObj,
						const LBM::BCData& BCInfo) {

		#pragma omp for collapse(2) schedule(static)
		for (size_t j = BCInfo.j_min; j <= BCInfo.j_max; ++j)  {
			for (size_t i = BCInfo.i_min; i <= BCInfo.i_max; ++i) {

				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid 
				size_t base = 9*n; 

				f[base + BCInfo.f_ref[0]] = f[base + BCInfo.f_inc[0]];
				f[base + BCInfo.f_ref[1]] = f[base + BCInfo.f_inc[1]];
				f[base + BCInfo.f_ref[2]] = f[base + BCInfo.f_inc[2]];
			}
		}
	}

	// Wall-tangent BC
	void tangentVelocityD2Q9(vector<double>& f,
							 const double* w,
							 const double cs2,
							 const LBM::CartesianGrid2D& gridObj,
							 const LBM::BCData& BCInfo) {

		#pragma omp for collapse(2) schedule(static)
		for (size_t j = BCInfo.j_min; j <= BCInfo.j_max; ++j)  {
			for (size_t i = BCInfo.i_min; i <= BCInfo.i_max; ++i) {

				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid 
				size_t base = 9*n; 

				// Calculate density at the current grid location
				double rho_ij = 0.0;
				for (size_t k = 0; k < 9; ++k) {
					double f_curr = f[9*n + k];
					rho_ij += f_curr;
				}

				f[base + BCInfo.f_ref[0]] = f[base + BCInfo.f_inc[0]] - 2*w[BCInfo.f_ref[0]]*rho_ij*BCInfo.U_wall/cs2;
				f[base + BCInfo.f_ref[1]] = f[base + BCInfo.f_inc[1]];
				f[base + BCInfo.f_ref[2]] = f[base + BCInfo.f_inc[2]] + 2*w[BCInfo.f_ref[2]]*rho_ij*BCInfo.U_wall/cs2;
			}
		}
	}


	// Legacy function
	void enforceD2Q9BCs(vector<double>& f,
						const double* w,
						const LBM::CartesianGrid2D& gridObj,
						const double U_lid,
						const double cs2,
						const double iyT) {

		// Bounce back on top
		#pragma omp for
		for (size_t i = 0; i < gridObj.Nx; ++i) {

			size_t n = i; // Current position on the flattened grid (j = 0)

			f[9*n + 2] = f[9*n + 6];
			f[9*n + 3] = f[9*n + 7];
			f[9*n + 4] = f[9*n + 8];
		}


		// Moving lid on top wall
		#pragma omp for
		for (size_t i = 0; i < gridObj.Nx; ++i) {
			size_t n = i + gridObj.Nx*iyT; // Current position on the flattened grid

			// Calculate density at current grid location
			double rho_ij = 0.0;
			for (size_t k = 0; k < 9; ++k) {
				double f_curr = f[9*n + k];
				rho_ij += f_curr;

			}

			// Update distribution function
			f[9*n + 6] = f[9*n + 2] - 2*w[6]*rho_ij*U_lid/cs2;
			f[9*n + 7] = f[9*n + 3];
			f[9*n + 8] = f[9*n + 4] + 2*w[8]*rho_ij*U_lid/cs2;
		}
	}
}
