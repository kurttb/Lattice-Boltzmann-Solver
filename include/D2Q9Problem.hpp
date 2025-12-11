#ifndef D2Q9_PROBLEM
#define D2Q9_PROBLEM

#include <vector>
#include <string>
#include "GridTypes.hpp"
#include "BoundaryConditions.hpp"
#include "GeneralProblem.hpp"


namespace LBM {

	class D2Q9Problem : public GeneralProblem
	{

		protected:

			// Computational grid structure
			CartesianGrid2D _gridObj; 

			// Boundary Conditions
			BCData _BCTop;
			BCData _BCBottom;
			BCData _BCRight;
			BCData _BCLeft;


		public:

			// Constructor/Destructor
			D2Q9Problem(const int Nx, const int Ny);
			~D2Q9Problem() override;

			// Pre-processing functions
			void setBC(const std::string& BCName, const std::string& BCType, const float uT = 0.0);
			// Run the simulation
			void runSimulation() override;

			// Post-process
			void writeOutput(std::string filePath) override;

			// Getters for testing
			std::vector<float> getRho() const;
			std::vector<float> getUx() const;
			std::vector<float> getUy() const;

		};

}



#endif
