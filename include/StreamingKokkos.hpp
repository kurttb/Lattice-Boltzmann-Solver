#ifndef STREAMING_KOKKOS
#define STREAMING_KOKKOS

#include <Kokkos_Core.hpp>

using vec2 = Kokkos::View<float**>;

struct ComputeStreaming {
	vec2_const f;
	vec2 fstream;
	const int Nx;
	const int Ny;

	ComputeStreaming(vec2_const f_, vec2 fstream_, const int Nx_, const int Ny_) :
		f(f_),
		fstream(fstream_),
		Nx(Nx_),
		Ny(Ny_) {}

	KOKKOS_INLINE_FUNCTION
	void operator() (const int n) const {
		static const int ex_stream[9] = {0, 1, 1, 0, -1, -1, -1,  0,  1};
		static const int ey_stream[9] = {0, 0, 1, 1,  1,  0, -1, -1, -1};

		int i = n % Nx;
		int j = n / Nx;

		for (int k = 0; k < 9; ++k) {
			int i_dest = (i + ex_stream[k] + Nx) % Nx;
			int j_dest = (j + ey_stream[k] + Ny) % Ny;

			int n_new = i_dest + Nx*j_dest;

			fstream(n_new, k) = f(n, k);
		}
	}
};


#endif
