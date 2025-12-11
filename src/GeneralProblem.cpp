#include <Kokkos_Core.hpp>
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "GeneralProblem.hpp"
#include "EquilibriumKokkos.hpp"
#include "ComputeStateKokkos.hpp"
#include "CollisionKokkos.hpp"
#include "StreamingKokkos.hpp"
#include "BCsKokkos.hpp"

namespace LBM {

	/////////// Preprocessor Functions
	
	// Constructor
	GeneralProblem::GeneralProblem() {
	}

	// Set Viscosity
	void GeneralProblem::setViscosity(const float nu) {
		if (nu < 0) {
			throw std::invalid_argument("Viscosity must be non-negative");
		}
		_nu = nu;
	}


	// Set initial condition
	void GeneralProblem::setIC(const float rho0, const float ux0, const float uy0) {
		if (rho0 <= 0) {
			throw std::invalid_argument("Density must be positive");
		}
		_rho0 = rho0;
		_ux0 = ux0;
		_uy0 = uy0;

	}



	// Set number of time steps
	void GeneralProblem::setNumTimeSteps(const int Nt) {
		_Nt = Nt;
	}

	// Set forces
	void GeneralProblem::setForces(const float Fx, const float Fy) {
		_Fx = Fx;
		_Fy = Fy;
	}


	/////////// Post-processing functions

	// Destructor
	GeneralProblem::~GeneralProblem() = default;
}

