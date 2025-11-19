#ifndef GRID_TYPES
#define GRID_TYPES

namespace LBM {

struct CartesianGrid2D {
	int Nx;
	int Ny;

	inline int idx(int i, int j) const {

		return i + j*Nx;
	}

};

}

#endif
