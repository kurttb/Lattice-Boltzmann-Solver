#ifndef VTK_WRITER
#define VTK_WRITER

#include <vector>
#include <string>

namespace LBM {
	void WriteVtk(const std::vector<double>& rho, 
				   const std::vector<double>& ux, 
				   const std::vector<double>& uy, 
				   const size_t Nx, const size_t Ny, 
				   const std::string& file_path, 
				   const std::string& pv_title);

}

#endif
