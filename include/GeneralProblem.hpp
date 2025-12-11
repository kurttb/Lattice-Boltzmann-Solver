#ifndef GENERAL_PROBLEM
#define GENERAL_PROBLEM

#include <vector>
#include <string>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"


namespace LBM {

	class GeneralProblem
	{

		protected:

			// Body Forces
			float _Fx = 0.0;
			float _Fy = 0.0;

			// Initial Conditions
			float _rho0 = 1.0;
			float _ux0 = 0.0;
			float _uy0 = 0.0;

			// Number of time steps
			int _Nt = 50000;

			// Fields
			std::vector<float> _rho;
			std::vector<float> _ux;
			std::vector<float> _uy;

			// Viscosity
			float _nu;

	public:

		// Constructor/Destructor
		GeneralProblem();
		virtual ~GeneralProblem();

		// Pre-processing functions
		void setIC(const float rho0, const float ux0, const float uy0);
		void setForces(const float Fx, const float Fy);
		void setNumTimeSteps(const int Nt);
		void setViscosity(const float nu);

		// Run the simulation
		virtual void runSimulation() = 0;

		// Post-process
		virtual void writeOutput(std::string filePath) = 0;

	};

}


#endif
