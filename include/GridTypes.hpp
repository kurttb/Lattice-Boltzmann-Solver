#ifndef GRID_TYPES
#define GRID_TYPES

namespace LBM {

struct CartesianGrid2D {
	size_t Nx;
	size_t Ny;

	int idx(size_t i, size_t j) const {

		return i + j*Nx;
	}

};

}

#endif
