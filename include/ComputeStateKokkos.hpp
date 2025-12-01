#ifndef COMPUTESTATEKOKKOS
#define COMPUTESTATEKOKKOS

#include <Kokkos_Core.hpp>

using vec1 = Kokkos::View<float*>;
using vec2 = Kokkos::View<float**>;
using vec2_const = Kokkos::View<const float**>;

struct ComputeState {
	vec1 rho;
	vec1 ux;
	vec1 uy;
	vec2_const f;
	const int* ex;
	const int* ey;
	const float Fx;
	const float Fy;
	const float tau;

	ComputeState(vec1 rho_, vec1 ux_, vec1 uy_, vec2_const f_, const int* ex_, const int* ey_, const float _Fx, const float _Fy, const float tau_) :
		rho(rho_),
		ux(ux_),
		uy(uy_),
		f(f_),
		ex(ex_),
		ey(ey_),
		Fx(_Fx),
		Fy(_Fy),
		tau(tau_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int n) const {
		float rho_ij = 0;
		float ux_ij = 0;
		float uy_ij = 0;

		for (int k = 0; k < 9; ++k) {
			float f_curr = f(n, k);
			rho_ij += f_curr;
			ux_ij += f_curr * ex[k];
			uy_ij += f_curr * ey[k];
		}

		ux(n) = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
		uy(n) = (uy_ij / rho_ij) + (Fy*tau / rho_ij);
		rho(n) = rho_ij;
	}
};

#endif

