#ifndef BCs_KOKKOS
#define BCs_KOKKOS

#include <Kokkos_Core.hpp>

using vec2 = Kokkos::View<float**>;

struct BounceBack {
	vec2 f;
	const int* f_inc;
	const int* f_ref;
	const int Nx;

	BounceBack(vec2 f_, const int* f_inc_, const int* f_ref_, const int Nx_) :
		f(f_),
		f_inc(f_inc_),
		f_ref(f_ref_),
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
	const int* f_inc;
	const int* f_ref;
	const int Nx;

	BounceBack(vec2 f_, const int* f_inc_, const int* f_ref_, const int Nx_) :
		f(f_),
		f_inc(f_inc_),
		f_ref(f_ref_),
		Nx(Nx_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int j, const int i) const {
		int n = i + Nx*j;

		f(n, f_ref[0]) = f(n, f_inc[0]);
		f(n, f_ref[1]) = f(n, f_inc[1]);
		f(n, f_ref[2]) = f(n, f_inc[2]);
	}
};


#endif
