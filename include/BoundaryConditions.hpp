// Header for boundary condition functions
#ifndef BCFUNCTIONS
#define BCFUNCTIONS

#include <vector>
#include "GridTypes.hpp"

using namespace std;

namespace LBM {
	void enforceD2Q9BCs(vector<double>& f,
						const double* w,
						const LBM::CartesianGrid2D& gridObj,
						const double U_lid,
						const double cs2,
						const double iyT);

}
#endif

