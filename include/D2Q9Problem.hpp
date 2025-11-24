#ifndef D2Q9_PROBLEM
#define D2Q9_PROBLEM

#include <vector>
#include <string>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"


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

			// Number of time steps
			size_t _Nt = 50000;

			// Fields
			vector<double> _rho;
			vector<double> _ux;
			vector<double> _uy;

			// Viscosity
			double _nu;

			// Distribution function
			vector<double> _f;

			// Define Lattice - Start at rest, go east, and move counterclockwise
			static constexpr int _ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
			static constexpr int _ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
			static constexpr double _w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distribution

			// Write path
			string _filePath;

			// Boundary Conditions
			BCData _BCTop;
			BCData _BCBottom;
			BCData _BCRight;
			BCData _BCLeft;


	public:

		// Constructor/Destructor
		D2Q9Problem(const size_t Nx, const size_t Ny);
		~D2Q9Problem();

		// Pre-processing functions
		void setIC(const double rho0, const double ux0, const double uy0);
		void setForces(const double Fx, const double Fy);
		void setNumTimeSteps(const size_t Nt);
		void setBC(const string& BCName, const string& BCType, const double uT = 0.0);
		void setViscosity(const double nu);

		// Run the simulation
		void runSimulation();

		// Post-process
		void writeOutput(string filePath);

	};

}







#endif
