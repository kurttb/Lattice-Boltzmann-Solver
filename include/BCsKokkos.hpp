#ifndef BCs_KOKKOS
#define BCs_KOKKOS

#include <Kokkos_Core.hpp>

using vec2 = Kokkos::View<float**>;

struct BounceBack {
	vec2 f;
	const int f_inc[3];
	const int f_ref[3];
	const int Nx;

	BounceBack(vec2 f_, const int f_inc_[3], const int f_ref_[3], const int Nx_) :
		f(f_),
		f_inc{f_inc_[0], f_inc_[1], f_inc_[2]},
		f_ref{f_ref_[0], f_ref_[1], f_ref_[2]},
		Nx(Nx_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int j, const int i) const {
		int n = i + Nx*j;

		f(n, f_ref[0]) = f(n, f_inc[0]);
		f(n, f_ref[1]) = f(n, f_inc[1]);
		f(n, f_ref[2]) = f(n, f_inc[2]);
	}
};




	struct TangentVelocity {
		vec2 f;
		const int f_inc[3];
		const int f_ref[3];
		const float U_wall;
		const float cs2;
		const int Nx;

		TangentVelocity(vec2 f_, const int f_inc_[3], const int f_ref_[3], const float U_wall_, const float cs2_, const int Nx_) :
			f(f_),
			f_inc{f_inc_[0], f_inc_[1], f_inc_[2]},
			f_ref{f_ref_[0], f_ref_[1], f_ref_[2]},
			U_wall(U_wall_),
			cs2(cs2_),
			Nx(Nx_) {}

		KOKKOS_INLINE_FUNCTION
		void operator() (const int j, const int i) const {
			static const float w_tan[9]  = {4.0f/9.0f,
									1.0f/9.0f, 1.0f/36.0f,
									1.0f/9.0f, 1.0f/36.0f,
									1.0f/9.0f, 1.0f/36.0f,
									1.0f/9.0f, 1.0f/36.0f};

			int n = i + Nx*j;

			// Calculate density at the current grid location
			float rho_ij = 0.0f;
			for (int k = 0; k < 9; ++k) {
				float f_curr = f(n, k);
				rho_ij += f_curr;
			}

			// Update distribution function
			f(n, f_ref[0]) = f(n, f_inc[0]) - 2*w_tan[f_ref[0]]*rho_ij*U_wall/cs2;
			f(n, f_ref[1]) = f(n, f_inc[1]);
			f(n, f_ref[2]) = f(n, f_inc[2]) + 2*w_tan[f_ref[2]]*rho_ij*U_wall/cs2;
		}
	};


#endif
