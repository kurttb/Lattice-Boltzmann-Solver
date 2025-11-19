// Header for Routines for recovering the state given the velocity distribution funciion

#ifndef RECONSTRUCT_STATE
#define RECONSTRUCT_STATE

#include <vector>
#include "GridTypes.hpp"

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
							  const double tau);

}
#endif


