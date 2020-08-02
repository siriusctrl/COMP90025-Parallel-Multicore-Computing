cmake_minimum_required(VERSION 3.16)
project(openmp_exercise)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Os")

set(CMAKE_CXX_STANDARD 14)
add_executable(main openmp.cpp openmp.h)
find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif()
