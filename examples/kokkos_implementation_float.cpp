// kokkos_implementation_float.cpp
// Port of naive_parallel.cpp to Kokkos (Float version)
// 11/27/25

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <Kokkos_Core.hpp>
#include "vtk_writer.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    Kokkos::initialize(argc, argv);
    {
        Kokkos::print_configuration(std::cout);

        // Define knobs
        // string flow_type = "Couette"; // Not used in logic directly except for Fx
        bool is_channel = false; // "Couette" by default
        
        constexpr size_t Nx = 1000; // Number of x coordinates
        constexpr size_t Ny = 1000; // Number of y coordinates
        constexpr size_t N = Nx*Ny; // Total number of grid nodes
        constexpr float Ma = 0.1f; // Mach number
        constexpr float Re = 100.0f; // Reynolds number
        constexpr size_t max_it = 50000; // Maximum number of iterations
        // constexpr float tol = 1e-4f; // Steady-state tolerance (unused in main loop)

        // Initial Condition
        constexpr float rho_init = 1.0f; // Initial density field
        constexpr float u_init[2] = {0.0f, 0.0f}; // Initial Velocity

        // Derive characteristics of the flow physics
        constexpr size_t L = Ny - 1; // Length of the domain in the lattice
        const float cs = 1.0f / sqrt(3.0f); // Speed of sound
        constexpr float cs2 = 1.0f / 3.0f; // Speed of sound squared
        const float U_lid = Ma*cs; // Lid velocity
        const float nu = (U_lid*L) / Re; // Kinematic viscosity
        const float tau = ( nu/(cs2) ) + 0.5f; // Relaxation parameter
        const float omega = 1.0f / tau;

        float Fx = 0;
        if (is_channel) {
            Fx += 1e-4f;
        }

        // Define Lattice constants on Device
        Kokkos::View<int[9]> ex("ex");
        Kokkos::View<int[9]> ey("ey");
        Kokkos::View<float[9]> w("w");

        // Initialize constants on Host and copy to Device
        {
            auto h_ex = Kokkos::create_mirror_view(ex);
            auto h_ey = Kokkos::create_mirror_view(ey);
            auto h_w = Kokkos::create_mirror_view(w);

            int ex_host[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
            int ey_host[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
            float w_host[9] = {4.0f/9.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f, 1.0f/9.0f, 1.0f/36.0f};

            for(int i=0; i<9; ++i) {
                h_ex(i) = ex_host[i];
                h_ey(i) = ey_host[i];
                h_w(i) = w_host[i];
            }
            Kokkos::deep_copy(ex, h_ex);
            Kokkos::deep_copy(ey, h_ey);
            Kokkos::deep_copy(w, h_w);
        }

        // Allocate Macroscopic Fields
        Kokkos::View<float*> rho("rho", N);
        Kokkos::View<float*> ux("ux", N);
        Kokkos::View<float*> uy("uy", N);

        // Initialize Macroscopic Fields
        Kokkos::parallel_for("InitMacro", N, KOKKOS_LAMBDA(const int n) {
            rho(n) = rho_init;
            ux(n) = u_init[0];
            uy(n) = u_init[1];
        });

        // Allocate Distribution Functions
        Kokkos::View<float**> f("f", N, 9);
        Kokkos::View<float**> fstream("fstream", N, 9);

        // Initialize to equilibrium
        Kokkos::parallel_for("InitEq", N, KOKKOS_LAMBDA(const int n) {
            float uxn = ux(n);
            float uyn = uy(n);
            float rho_n = rho(n);
            float u_sq_ind = uxn*uxn + uyn*uyn;

            for(int k = 0; k < 9; ++k) {
                float e_dot_u = uxn*static_cast<float>(ex(k)) + uyn*static_cast<float>(ey(k));
                f(n, k) = w(k) * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
            }
        });

        // Boundary Indices
        size_t iyT = Ny - 1;

        Kokkos::Timer timer;

        // Start Update Loop
        for (size_t it = 0; it < max_it; ++it) {

            // Compute macroscopic quantities from distribution
            Kokkos::parallel_for("ComputeMacro", N, KOKKOS_LAMBDA(const int n) {
                float rho_ij = 0;
                float ux_ij = 0;
                float uy_ij = 0;

                for (int k = 0; k < 9; ++k) {
                    float f_curr = f(n, k);
                    rho_ij += f_curr;
                    ux_ij += f_curr * ex(k);
                    uy_ij += f_curr * ey(k);
                }

                ux(n) = (ux_ij / rho_ij) + (Fx*tau / rho_ij);
                uy(n) = uy_ij / rho_ij;
                rho(n) = rho_ij;
            });

            // Collision Step
            Kokkos::parallel_for("Collision", N, KOKKOS_LAMBDA(const int n) {
                float uxn = ux(n);
                float uyn = uy(n);
                float rho_n = rho(n);
                float u_sq_ind = uxn*uxn + uyn*uyn;

                for(int k = 0; k < 9; ++k) {
                    float e_dot_u = uxn*static_cast<float>(ex(k)) + uyn*static_cast<float>(ey(k));
                    float f_curr = f(n, k);
                    float feq_curr = w(k) * rho_n * (1.0f + 3.0f*e_dot_u + 4.5f*e_dot_u*e_dot_u - 1.5f*u_sq_ind);
                    f(n, k) = f_curr - omega * (f_curr - feq_curr);
                }
            });

            // Advection step
            Kokkos::parallel_for("Streaming", N, KOKKOS_LAMBDA(const int n) {
                int i = n % Nx;
                int j = n / Nx;

                for(int k = 0; k < 9; ++k) {
                    int i_dest = (i + ex(k) + Nx) % Nx;
                    int j_dest = (j + ey(k) + Ny) % Ny;
                    
                    int n_new = i_dest + Nx*j_dest;
                    
                    fstream(n_new, k) = f(n, k);
                }
            });

            // Swap f and fstream
            std::swap(f, fstream);

            // Bounce-back on bottom wall (j=0)
            Kokkos::parallel_for("BounceBackBottom", Nx, KOKKOS_LAMBDA(const int i) {
                int n = i; // j=0
                f(n, 2) = f(n, 6);
                f(n, 3) = f(n, 7);
                f(n, 4) = f(n, 8);
            });

            // Moving lid on top wall (j=Ny-1)
            Kokkos::parallel_for("MovingLidTop", Nx, KOKKOS_LAMBDA(const int i) {
                int n = i + Nx*iyT;
                
                float rho_ij = 0.0f;
                for (int k = 0; k < 9; ++k) {
                    rho_ij += f(n, k);
                }

                f(n, 6) = f(n, 2) - 2.0f*w(6)*rho_ij*U_lid/cs2;
                f(n, 7) = f(n, 3);
                f(n, 8) = f(n, 4) + 2.0f*w(8)*rho_ij*U_lid/cs2;
            });
        }

        Kokkos::fence();
        double time = timer.seconds();
        printf("Runtime: %f s\n", time);

        // Write output
        // Copy back to host std::vector for the write_vtk function
        // NOTE: write_vtk expects double, so we keep these as double and cast
        std::vector<double> h_rho_vec(N);
        std::vector<double> h_ux_vec(N);
        std::vector<double> h_uy_vec(N);

        auto h_rho = Kokkos::create_mirror_view(rho);
        auto h_ux = Kokkos::create_mirror_view(ux);
        auto h_uy = Kokkos::create_mirror_view(uy);

        Kokkos::deep_copy(h_rho, rho);
        Kokkos::deep_copy(h_ux, ux);
        Kokkos::deep_copy(h_uy, uy);

        for(size_t i=0; i<N; ++i) {
            h_rho_vec[i] = static_cast<double>(h_rho(i));
            h_ux_vec[i] = static_cast<double>(h_ux(i));
            h_uy_vec[i] = static_cast<double>(h_uy(i));
        }

        std::string file_path = "./tmp/paraview_kokkos_float.vtk";
        std::string pv_title = "LBM Field Float";
        write_vtk(h_rho_vec, h_ux_vec, h_uy_vec, Nx, Ny, file_path, pv_title);
    }
    Kokkos::finalize();

    return 0;
}
