// Header for boundary condition functions
#ifndef BCFUNCTIONS
#define BCFUNCTIONS

#include <vector>
#include "GridTypes.hpp"
#include <string>

using namespace std;

namespace LBM {

	// Structure for storing boundary condition data
	struct BCData {
		string BCType; // Name of boundary condition
		int f_inc[3]; // Incident distribution function IDs
		int f_ref[3]; // Reflected distribution function IDs
		double U_wall; // Wall velocity (if applicable)
		int i_min; // Minimum i for the boundary condtion
		int i_max; // Maximum i for the boundary condition
		int j_min; // Minimum j for the boundary condtion
		int j_max; // Maximum j for the boundary condition
	};

	// Bounce-Back BC
	void bounceBackD2Q9(vector<double>& f,
						const LBM::CartesianGrid2D& gridObj,
						const LBM::BCData& BCInfo);

	// Wall-tangent BC
	void tangentVelocityD2Q9(vector<double>& f,
							 const double* w,
							 const double cs2,
							 const LBM::CartesianGrid2D& gridObj,
							 const LBM::BCData& BCInfo);
	
	// Initial Hard-Coded BC
	void enforceD2Q9BCs(vector<double>& f,
						const double* w,
						const LBM::CartesianGrid2D& gridObj,
						const double U_lid,
						const double cs2,
						const double iyT);

}
#endif

