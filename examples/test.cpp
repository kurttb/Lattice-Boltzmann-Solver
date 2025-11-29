// config_test.cpp
#include <Kokkos_Core.hpp>
#include <iostream>

int main() {
    Kokkos::initialize();
    {
        Kokkos::print_configuration(std::cout);
    }
    Kokkos::finalize();
    return 0;
}

