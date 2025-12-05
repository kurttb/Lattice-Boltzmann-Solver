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
	const float Cs;

	ComputeCollision(vec1_const rho_, vec1_const ux_, vec1_const uy_, vec2 f_, const float omega_, const float Cs_) :
		rho(rho_),
		ux(ux_),
		uy(uy_),
		f(f_),
		omega(omega_),
		Cs(Cs_) {}

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

		float f_neq[9];
		float feq[9];

		// Compute Equilibrium and Non-Equilibrium parts
		for (int k = 0; k < 9; ++k) {
			float e_dot_u = uxn*static_cast<float>(ex_coll[k]) + uyn*static_cast<float>(ey_coll[k]);
			feq[k] = w_coll[k] * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
			f_neq[k] = f(n, k) - feq[k];
		}

		float omega_eff = omega;

		if (Cs > 0.0f) {
			// Compute Momentum Flux Tensor Q_alpha_beta
			float Qxx = 0.0f;
			float Qxy = 0.0f;
			float Qyy = 0.0f;

			for (int k = 0; k < 9; ++k) {
				float ex_k = static_cast<float>(ex_coll[k]);
				float ey_k = static_cast<float>(ey_coll[k]);
				
				Qxx += ex_k * ex_k * f_neq[k];
				Qxy += ex_k * ey_k * f_neq[k];
				Qyy += ey_k * ey_k * f_neq[k];
			}

			// Compute P = sqrt(2 * sum(Q_ab * Q_ab))
			// sum(Q_ab * Q_ab) = Qxx^2 + Qyy^2 + 2*Qxy^2
			float P = sqrt(2.0f * (Qxx*Qxx + Qyy*Qyy + 2.0f*Qxy*Qxy));

			// Compute effective relaxation time
			float tau0 = 1.0f / omega;
			float S_param = 18.0f * Cs * Cs * P / rho_n; // 18 * Cs^2 * P / rho
			// tau_eff = (tau0 + sqrt(tau0^2 + S_param)) / 2
			float tau_eff = 0.5f * (tau0 + sqrt(tau0*tau0 + S_param));
			
			omega_eff = 1.0f / tau_eff;
		}

		for (int k = 0; k < 9; ++k) {
			f(n, k) = f(n, k) - omega_eff * f_neq[k];
		}
	}
};


#endif
