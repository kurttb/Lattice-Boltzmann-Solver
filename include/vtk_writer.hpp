#include <vector>
#include <string>

#ifndef VTK_WRITER
#define VTK_WRITER

void write_vtk(const std::vector<double>& rho, const std::vector<double>& ux, const std::vector<double>& uy, const int Nx, const int Ny, const std::string& file_path, const std::string& pv_title);




#endif
