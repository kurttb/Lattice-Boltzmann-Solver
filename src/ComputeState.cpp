// Implementation for Routines for recovering the state given the velocity distribution funciion

#include <vector>
#include "GridTypes.hpp"
#include "ComputeState.hpp"

using namespace std;

namespace LBM {
	// Reconstruct D2Q9 state (uniform cartesian grid)
	void D2Q9ReconstructState(vector<double>& rho, 
							  vector<double>& ux, 
							  vector <double>& uy, 
							  const vector<double>& f, 
							  const int* ex,
							  const int* ey,
							  const LBM::CartesianGrid2D& gridObj,
							  const double Fx,
							  const double Fy,
							  const double tau) {

		#pragma omp for collapse(2) schedule(static)
		for(size_t j = 0; j < gridObj.Ny; ++j) {
			for (size_t i = 0; i < gridObj.Nx; ++i) {
				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
				double rho_ij = 0;
				double ux_ij = 0;
				double uy_ij = 0;

				for (size_t k = 0; k < 9; ++k) {
					double f_curr = f[9*n + k];
					rho_ij += f_curr;
					ux_ij += f_curr * ex[k];
					uy_ij += f_curr * ey[k];
				}

				ux[n] = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
				uy[n] = (uy_ij / rho_ij) + (Fy*tau / rho_ij);
				rho[n] = rho_ij;
			}
		}

	}
}






