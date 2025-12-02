#ifndef COLLISION_KOKKOS
#define COLLISION_KOKKOS

#include <Kokkos_Core.hpp>

using vec1_const = Kokkos::View<const float*>;
using vec2 = Kokkos::View<float**>;

struct ComputeCollision {
	vec1_const rho;
	vec1_const ux;
	vec1_const uy;
	vec2 f;
	const float omega;

	ComputeCollision(vec1_const rho_, vec1_const ux_, vec1_const uy_, vec2 f_, const float omega_) :
		rho(rho_),
		ux(ux_),
		uy(uy_),
		f(f_),
		omega(omega_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int n) const {
		static const int   ex_coll[9] = {0, 1, 1, 0, -1, -1, -1,  0,  1};
		static const int   ey_coll[9] = {0, 0, 1, 1,  1,  0, -1, -1, -1};
		static const float w_coll[9]  = {4.0f/9.0f,
                                            1.0f/9.0f, 1.0f/36.0f,
                                            1.0f/9.0f, 1.0f/36.0f,
                                            1.0f/9.0f, 1.0f/36.0f,
                                            1.0f/9.0f, 1.0f/36.0f};
		float uxn = ux(n);
		float uyn = uy(n);
		float rho_n = rho(n);
		float u_sq_ind = uxn*uxn + uyn*uyn;

		for (int k = 0; k < 9; ++k) {
			float e_dot_u = uxn*static_cast<float>(ex_coll[k]) + uyn*static_cast<float>(ey_coll[k]);
			float f_curr = f(n, k);
			float feq_curr = w_coll[k] * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
			f(n, k) = f_curr - omega * (f_curr - feq_curr);
		}
	}
};


#endif
