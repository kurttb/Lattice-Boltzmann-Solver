// Header file for different collision models
#ifndef COLLISIONS
#define COLLISIONS

#include <vector>
#include "GridTypes.hpp"

using namespace std;

namespace LBM{

	// Bhatnagar-Gross-Krook Collision model for the D2Q9 lattice
	void D2Q9BGKCollision(const vector<double>& rho,
					   const vector<double>& ux,
					   const vector<double>& uy,
					   vector<double>& f,
					   const int* ex,
					   const int* ey,
					   const double* w,
					   const LBM::CartesianGrid2D& gridObj,
					   const double omega);

}
#endif
