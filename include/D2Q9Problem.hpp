#ifndef D2Q9_PROBLEM
#define D2Q9_PROBLEM

#include <vector>
#include <string>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"


namespace LBM {

	class D2Q9Problem
	{

		private:

			// Computational grid structure
			CartesianGrid2D _gridObj; 

			// Body Forces
			float _Fx = 0.0;
			float _Fy = 0.0;

			// Initial Conditions
			float _rho0 = 1.0;
			float _ux0 = 0.0;
			float _uy0 = 0.0;

			// Number of time steps
			size_t _Nt = 50000;

			// Fields
			vector<float> _rho;
			vector<float> _ux;
			vector<float> _uy;

			// Viscosity
			float _nu;

			// Distribution function
			vector<float> _f;

			// Define Lattice - Start at rest, go east, and move counterclockwise
			static constexpr int _ex[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
			static constexpr int _ey[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
			static constexpr float _w[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0}; // Weights for Maxwellian Distribution

			// Write path
            std::string _filePath;

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
		void setIC(const float rho0, const float ux0, const float uy0);
		void setForces(const float Fx, const float Fy);
		void setNumTimeSteps(const size_t Nt);
		void setBC(const std::string& BCName, const std::string& BCType, const float uT = 0.0);
		void setViscosity(const float nu);

		// Run the simulation
		void runSimulation();

		// Post-process
		void writeOutput(std::string filePath);

	};

}







#endif
