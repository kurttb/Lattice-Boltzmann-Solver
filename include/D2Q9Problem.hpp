#ifndef D2Q9_PROBLEM
#define D2Q9_PROBLEM

#include <vector>
#include <string>
#include "GridTypes.hpp"


using namespace std;

namespace LBM {

	class D2Q9Problem
	{

		private:

			// Computational grid structure
			CartesianGrid2D _gridObj; 

			// Body Forces
			double _Fx = 0.0;
			double _Fy = 0.0;

			// Initial Conditions (scalars)
			double _rho0 = 1.0;
			double _ux0 = 0.0;
			double _uy0 = 0.0;

			// Number of time steps
			int _Nt = 50000;

			// Fields
			vector<double> _rho;
			vector<double> _ux;
			vector<double> _uy;

			// Distribution function
			vector<double> _f;

			// Define Lattice - Start at rest, go east, and move counterclockwise
			static constexpr int _ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
			static constexpr int _ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
			static constexpr double _w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distribution

			// Write path
			string _filePath;

			// Boundary Conditions
			struct _BCData {

				// BC Types
				string _TBCType; // Top BC type
				string _BBCType; // Bottom BC type
				string _RBCType; // Right BC type
				string _LBCType; // Left BC type

				// BC values
				double uxT; // Top ux
				double uyT; // Top uy
				double uxB; // Bottom ux
				double uyB; // Bottom uy
			};


	public:

		// Constructor/Destructor
		D2Q9Problem(const int Nx, const int Ny);
		~D2Q9Problem();

		// Pre-processing functions
		void setIC(const double rho0, const double ux0, const double uy0);
		void setForces(const double Fx, const double Fy);
		void setBC(const string& BCLabel, const string& BCType, const double ux = 0.0, const double uy = 0.0);
		void setNumTimeSteps(const int Nt);

		// Run the simulation
		void runSimulation();

		// Post-process
		void writeOutput(string filePath);

	};

}







#endif
