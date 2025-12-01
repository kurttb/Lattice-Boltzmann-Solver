#ifndef STREAMING_KOKKOS
#define STREAMING_KOKKOS

#include <Kokkos_Core.hpp>

using vec2 = Kokkos::View<float**>;
using vec2_const = Kokkos::View<const float**>;

struct ComputeStreaming {
	vec2_const f;
	vec2 fstream;
	const int* ex;
	const int* ey;
	const int Nx;
	const int Ny;

	ComputeStreaming(vec2_const f_, vec2 fstream_, const int* ex_, const int* ey_, const int Nx_, const int Ny_) :
		f(f_),
		fstream(fstream_),
		ex(ex_),
		ey(ey_),
		Nx(Nx_),
		Ny(Ny_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int n) const {
		int i = n % Nx;
		int j = n / Nx;

		for (int k = 0; k < 9; ++k) {
			int i_dest = (i + ex[k] + Nx) % Nx;
			int j_dest = (j + ey[k] + Ny) % Ny;

			int n_new = i_dest + Nx*j_dest;

			fstream(n_new, k) = f(n, k);
		}
	}
};


#endif
