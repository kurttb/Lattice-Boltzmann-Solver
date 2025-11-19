// Implementation of Streaming step functions
#include <vector>
#include "GridTypes.hpp"
#include "Streaming.hpp"

using namespace std;

namespace LBM {
	void D2Q9Stream(const vector<double>& f, 
					vector<double>& fstream, 
					const int* ex,
					const int* ey,
					LBM::CartesianGrid2D& gridObj) {

		#pragma omp for collapse(2) schedule(static)
		for(size_t j = 1; j < gridObj.Ny - 1; ++j) {
			for(size_t i = 1; i < gridObj.Nx - 1; ++i) {
				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
				size_t base = 9*n;

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = i + ex[k];
					size_t j_new = j + ey[k];
					size_t n_new = i_new + gridObj.Nx*j_new; // Location of point to be advected to

					fstream[9*n_new + k] = f[base + k];

				}
			}
		}


		// Streaming left/right
		#pragma omp for
		for(size_t j = 1; j < gridObj.Ny - 1; ++j) {
			for(auto i: { size_t(0), size_t(gridObj.Nx - 1) }) {
				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
				size_t base = 9*n;

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = (i + ex[k] + gridObj.Nx) % gridObj.Nx;
					size_t j_new = j + ey[k];
					size_t n_new = i_new + gridObj.Nx*j_new; // Location of point to be advected to

					fstream[9*n_new + k] = f[base + k];

				}
			}
		}


		// Streaming top/bottom
		#pragma omp for
		for(size_t i = 1; i < gridObj.Nx - 1; ++i) {
			for(auto j: { size_t(0), size_t(gridObj.Ny - 1) }) {
				size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
				size_t base = 9*n;

				for(size_t k = 0; k < 9; ++k) {
					size_t i_new = i + ex[k];
					size_t j_new = (j + ey[k] + gridObj.Ny) % gridObj.Ny;

					size_t n_new = i_new + gridObj.Nx*j_new; // Location of point to be advected to

					fstream[9*n_new + k] = f[base + k];

				}
			}
		}



		// Stream the corners
		#pragma omp single
		{
			for(auto i: {size_t(0), size_t(gridObj.Nx - 1)}) {
				for(auto j: { size_t(0), size_t(gridObj.Ny - 1) }) {
					size_t n = i + gridObj.Nx*j; // Current position on the flattened grid
					size_t base = 9*n;

					for(size_t k = 0; k < 9; ++k) {
						size_t i_new = (i + ex[k] + gridObj.Nx) % gridObj.Nx;
						size_t j_new = (j + ey[k] + gridObj.Ny) % gridObj.Ny;
						size_t n_new = i_new + gridObj.Nx*j_new; // Location of point to be advected to

						fstream[9*n_new + k] = f[base + k];

					}
				}
			}

		}

	}
}
