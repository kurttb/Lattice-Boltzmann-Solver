// Header for boundary condition functions
#ifndef BCFUNCTIONS
#define BCFUNCTIONS

#include "GridTypes.hpp"
#include <string>

using namespace std;

namespace LBM {

	// Structure for storing boundary condition data
	struct BCData {
		string BCType; // Name of boundary condition
		int f_inc[3]; // Incident distribution function IDs
		int f_ref[3]; // Reflected distribution function IDs
		float U_wall; // Wall velocity (if applicable)
		int i_min; // Minimum i for the boundary condtion
		int i_max; // Maximum i for the boundary condition
		int j_min; // Minimum j for the boundary condtion
		int j_max; // Maximum j for the boundary condition
	};
}
#endif

