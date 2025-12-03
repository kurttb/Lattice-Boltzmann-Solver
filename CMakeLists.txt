cmake_minimum_required(VERSION 3.15)

project(LBM VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define Variables for Kokkos Backend
set(Kokkos_ROOT "/usr" CACHE PATH "Root directory of Kokkos installation")

# Find Kokkos
find_package(Kokkos REQUIRED)


# Set Location for Libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)


# Add Library and its source files
add_library(LBM STATIC
		src/D2Q9Problem.cpp
		src/VtkWriter.cpp
)

# Link LBM library to kokkos
target_link_libraries(LBM
	PRIVATE
		Kokkos::kokkos
)

# Add include directory
target_include_directories(LBM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set location for exectuables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# Add exectuables
file(GLOB files "examples/*_case.cpp")
foreach(file ${files})
	
	get_filename_component(name ${file} NAME_WE) # Get name without extension
	add_executable(${name}) # Make an executable

	target_sources(${name}
		PRIVATE
			${file}
	) # Add source file for target

	target_link_libraries(${name}
		PRIVATE
			LBM
	) # Link against the LBM library


	# -I flag for the specific target
	target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

endforeach()



# Add compiler flag macro
macro(add_compiler_flag FLAG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
endmacro()


# Add compiler flags
add_compiler_flag("-O3")
add_compiler_flag("-ffast-math")
add_compiler_flag("-Wall")
add_compiler_flag("-Wextra")
add_compiler_flag("-Wpedantic")
