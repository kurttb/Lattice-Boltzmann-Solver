// kokkos_implementation.cpp
// Port of naive_parallel.cpp to Kokkos
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
        constexpr double Ma = 0.1; // Mach number
        constexpr double Re = 100; // Reynolds number
        constexpr size_t max_it = 50000; // Maximum number of iterations
        // constexpr double tol = 1e-4; // Steady-state tolerance (unused in main loop)

        // Initial Condition
        constexpr double rho_init = 1.0; // Initial density field
        constexpr double u_init[2] = {0.0, 0.0}; // Initial Velocity

        // Derive characteristics of the flow physics
        constexpr size_t L = Ny - 1; // Length of the domain in the lattice
        const double cs = 1.0 / sqrt(3); // Speed of sound
        constexpr double cs2 = 1.0 / 3.0; // Speed of sound squared
        const double U_lid = Ma*cs; // Lid velocity
        const double nu = (U_lid*L) / Re; // Kinematic viscosity
        const double tau = ( nu/(cs2) ) + 0.5; // Relaxation parameter
        const double omega = 1.0 / tau;

        double Fx = 0;
        if (is_channel) {
            Fx += 1e-4;
        }

        // Define Lattice constants on Device
        Kokkos::View<int[9]> ex("ex");
        Kokkos::View<int[9]> ey("ey");
        Kokkos::View<double[9]> w("w");

        // Initialize constants on Host and copy to Device
        {
            auto h_ex = Kokkos::create_mirror_view(ex);
            auto h_ey = Kokkos::create_mirror_view(ey);
            auto h_w = Kokkos::create_mirror_view(w);

            int ex_host[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
            int ey_host[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
            double w_host[9] = {4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0};

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
        Kokkos::View<double*> rho("rho", N);
        Kokkos::View<double*> ux("ux", N);
        Kokkos::View<double*> uy("uy", N);

        // Initialize Macroscopic Fields
        Kokkos::parallel_for("InitMacro", N, KOKKOS_LAMBDA(const int n) {
            rho(n) = rho_init;
            ux(n) = u_init[0];
            uy(n) = u_init[1];
        });

        // Allocate Distribution Functions
        // Using LayoutRight (N, 9) to match original indexing logic [9*n + k] -> (n, k)
        // But for GPU coalescing, (n, k) with LayoutLeft (default on GPU) means stride 1 is on n.
        // If we use View<double**> f("f", N, 9), default layout on GPU is LayoutLeft.
        // So f(n, k) -> memory at n + k*N.
        // Original code: f[9*n + k] -> memory at 9*n + k. This is LayoutRight (row-major).
        // Let's use Kokkos default layout.
        Kokkos::View<double**> f("f", N, 9);
        Kokkos::View<double**> fstream("fstream", N, 9);

        // Initialize to equilibrium
        Kokkos::parallel_for("InitEq", N, KOKKOS_LAMBDA(const int n) {
            double uxn = ux(n);
            double uyn = uy(n);
            double rho_n = rho(n);
            double u_sq_ind = uxn*uxn + uyn*uyn;

            for(int k = 0; k < 9; ++k) {
                double e_dot_u = uxn*static_cast<double>(ex(k)) + uyn*static_cast<double>(ey(k));
                f(n, k) = w(k) * rho_n * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
            }
        });

        // Boundary Indices
        // size_t ixL = 0;
        // size_t ixR = Nx - 1;
        // size_t iyB = 0;
        size_t iyT = Ny - 1;

        Kokkos::Timer timer;

        // Start Update Loop
        for (size_t it = 0; it < max_it; ++it) {

            // Compute macroscopic quantities from distribution
            Kokkos::parallel_for("ComputeMacro", N, KOKKOS_LAMBDA(const int n) {
                double rho_ij = 0;
                double ux_ij = 0;
                double uy_ij = 0;

                for (int k = 0; k < 9; ++k) {
                    double f_curr = f(n, k);
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
                double uxn = ux(n);
                double uyn = uy(n);
                double rho_n = rho(n);
                double u_sq_ind = uxn*uxn + uyn*uyn;

                for(int k = 0; k < 9; ++k) {
                    double e_dot_u = uxn*static_cast<double>(ex(k)) + uyn*static_cast<double>(ey(k));
                    double f_curr = f(n, k);
                    double feq_curr = w(k) * rho_n * (1 + 3*e_dot_u + 4.5*e_dot_u*e_dot_u - 1.5*u_sq_ind);
                    f(n, k) = f_curr - omega * (f_curr - feq_curr);
                }
            });

            // Advection step
            // We can use a single kernel for all streaming if we handle boundaries carefully.
            // Or we can replicate the structure of the original code.
            // Original code splits: Interior, Left/Right, Top/Bottom, Corners.
            // Let's try to do it in one go with conditionals, or use MDRangePolicy.
            // Using conditionals inside one kernel is easier to write and often fine.
            
            Kokkos::parallel_for("Streaming", N, KOKKOS_LAMBDA(const int n) {
                int i = n % Nx;
                int j = n / Nx;

                for(int k = 0; k < 9; ++k) {
                    // Determine next position
                    int i_new = i + ex(k);
                    int j_new = j + ey(k);

                    // Periodic Boundary Conditions (wrap around)
                    // The original code does explicit periodic wrapping for boundaries.
                    // Interior: i_new = i + ex[k], j_new = j + ey[k]
                    // Left/Right: i_new = (i + ex[k] + Nx) % Nx
                    // Top/Bottom: j_new = (j + ey[k] + Ny) % Ny
                    
                    // We can just apply modulo arithmetic everywhere for periodic boundaries
                    // But wait, the original code has "Advection step on domain (Excluding top and right boundaries)"?
                    // No, "Excluding top and right boundaries" comment in original code seems to refer to the loop limits `j < Ny - 1`, `i < Nx - 1`.
                    // But then it has "Streaming left/right" and "Streaming top/bottom" and "Stream the corners".
                    // It seems it implements full periodic boundaries by parts.
                    
                    int i_dest = (i + ex(k) + Nx) % Nx;
                    int j_dest = (j + ey(k) + Ny) % Ny;
                    
                    int n_new = i_dest + Nx*j_dest;
                    
                    fstream(n_new, k) = f(n, k);
                }
            });

            // Swap f and fstream
            // We can't swap Views inside the device, but we can swap the handles on the host.
            // However, we need to make sure the next iteration uses the swapped views.
            // Since we are in a host loop, we can just swap the View objects.
            std::swap(f, fstream);

            // Bounce-back on bottom wall (j=0)
            Kokkos::parallel_for("BounceBackBottom", Nx, KOKKOS_LAMBDA(const int i) {
                int n = i; // j=0
                // f[9*n + 2] = f[9*n + 6];
                // f[9*n + 3] = f[9*n + 7];
                // f[9*n + 4] = f[9*n + 8];
                f(n, 2) = f(n, 6);
                f(n, 3) = f(n, 7);
                f(n, 4) = f(n, 8);
            });

            // Moving lid on top wall (j=Ny-1)
            Kokkos::parallel_for("MovingLidTop", Nx, KOKKOS_LAMBDA(const int i) {
                int n = i + Nx*iyT;
                
                double rho_ij = 0.0;
                for (int k = 0; k < 9; ++k) {
                    rho_ij += f(n, k);
                }

                f(n, 6) = f(n, 2) - 2*w(6)*rho_ij*U_lid/cs2;
                f(n, 7) = f(n, 3);
                f(n, 8) = f(n, 4) + 2*w(8)*rho_ij*U_lid/cs2;
            });
        }

        Kokkos::fence();
        double time = timer.seconds();
        printf("Runtime: %f s\n", time);

        // Write output
        // Copy back to host std::vector for the write_vtk function
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
            h_rho_vec[i] = h_rho(i);
            h_ux_vec[i] = h_ux(i);
            h_uy_vec[i] = h_uy(i);
        }

        std::string file_path = "./tmp/paraview_kokkos.vtk";
        std::string pv_title = "LBM Field";
        write_vtk(h_rho_vec, h_ux_vec, h_uy_vec, Nx, Ny, file_path, pv_title);
    }
    Kokkos::finalize();

    return 0;
}
