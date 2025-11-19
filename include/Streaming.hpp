// Header file for streaming functions
#ifndef STREAMING
#define STREAMING

#include <vector>
#include "GridTypes.hpp"
#include "Streaming.hpp"

using namespace std;


namespace LBM {

	void D2Q9Stream(const vector<double>& f,
					vector<double>& fstream,
					const int* ex,
					const int* ey,
					LBM::CartesianGrid2D& gridObj);
}

#endif
