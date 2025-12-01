// Header file for Kokkos equilibrium calculation
#ifndef CALC_EQ_KOKKOS
#define CALC_EQ_KOKKOS

#include <Kokkos_Core.hpp>

using vec1 = Kokkos::View<float*>;
using vec1_const = Kokkos::View<const float*>;
using vec2 = Kokkos::View<float**>;

struct CalcEq {
	vec1_const rho;
	vec1_const ux;
	vec1_const uy;
	vec2 f;
	const int* ex;
	const int* ey;
	const float* w;

	CalcEq(vec1_const rho_, vec1_const ux_, vec1_const uy_, vec2 f_, const int* ex_, const int* ey_, const float* w_) :
		rho(rho_),
		ux(ux_),
		uy(uy_),
		f(f_),
		ex(ex_),
		ey(ey_),
		w(w_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int n) const {
		float uxn = ux(n);
		float uyn = uy(n);
		float rho_n = rho(n);
		float u_sq_ind = uxn*uxn + uyn*uyn;

		for (int k = 0; k < 9; ++k) {
			float e_dot_u = uxn*static_cast<float>(ex[k]) + uyn*static_cast<float>(ey[k]);
			f(n, k) = w[k] * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
		}
	}
};





#endif 
