// Implementation of boundary condition functions
#include <vector>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"

using namespace std;

namespace LBM {
	void enforceD2Q9BCs(const vector<double>& rho, 
						vector<double>& f,
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
