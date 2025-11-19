// Implementation for different collision models
#include <vector>
#include "GridTypes.hpp"
#include "Collisions.hpp"

using namespace std;

namespace LBM {
	void D2Q9BGKCollision(const vector<double>& rho,
					   const vector<double>& ux,
					   const vector<double>& uy,
					   vector<double>& f,
					   const int* ex,
					   const int* ey,
					   const double* w,
					   const LBM::CartesianGrid2D& gridObj,
					   const double omega) {

		#pragma omp for collapse(2) schedule(static)
		for(size_t j = 0; j < gridObj.Ny; ++j) {
			for(size_t i = 0; i < gridObj.Nx; ++i) {
				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
				double uxn = ux[n]; // Extract x velocity component
				double uyn = uy[n]; // Extract y velocity component
				double u_sq_ind = ux[n]*ux[n] + uy[n]*uy[n]; // Square velocity at the current grid node
				size_t base = 9*n;

				for(size_t k = 0; k < 9; ++k) {
					double e_dot_u = uxn*static_cast<double>(ex[k]) + uyn*static_cast<double>(ey[k]);
					double f_curr = f[base + k];
					double feq_curr = w[k] * rho[n] * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind); // Calculate equilbrium distribution
					f[base + k] = f_curr - omega * (f_curr - feq_curr);
				}
			}
		}




	}

}
