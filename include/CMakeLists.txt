hello(${CMAKE_CURRENT_SOURCE_DIR})

# This is like adding -I <this directory> to compiler flags
set(LBM_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "Include directory for LBM headers")
