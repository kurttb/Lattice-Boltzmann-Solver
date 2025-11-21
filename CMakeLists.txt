cmake_minimum_required(VERSION 3.23)

project(LBM VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Set Location for Libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)


# Add Library and its source files
add_library(LBM SHARED
		src/BoundaryConditions.cpp
		src/Collisions.cpp
		src/ComputeState.cpp
		src/D2Q9Problem.cpp
		src/Streaming.cpp
		src/VtkWriter.cpp
)

# Set location for exectuables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Add exectuables (example)
add_executable(Couette)


# Add source files to executable
target_sources(Couette
	PRIVATE
		examples/Couette.cpp
)

# Link Exectuable with the library
target_link_libraries(Couette
	PRIVATE
		LBM
)

# Add include directory
target_include_directories(LBM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(Couette PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


# Add compiler flag macro
macro(add_compiler_flag FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
endmacro()


# Add compiler flags
add_compiler_flag("-O3")
add_compiler_flag("-ffast-math")
add_compiler_flag("-march=native")
add_compiler_flag("-Wall")
add_compiler_flag("-Wextra")
add_compiler_flag("-Wpedantic")


# Find OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(LBM PRIVATE OpenMP::OpenMP_CXX)



# Write Directory
#set_target_properties(LBM 
#	PROPERTIES
#		LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
#)
