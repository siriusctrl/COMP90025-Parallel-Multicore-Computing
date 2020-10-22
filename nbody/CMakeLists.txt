cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_COMPILER "mpicxx")

project(nbody_nlogn)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

add_executable(nlogn nlogn.cpp)

# link MPI
target_link_libraries(nlogn ${MPI_LIBRARIES})

# add openmp dependency
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nlogn PUBLIC OpenMP::OpenMP_CXX)
endif()

set(MPI_COMPILE_FLAGS "-O3")