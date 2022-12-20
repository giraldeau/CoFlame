# We require CMake 3.10 for the newer FindMPI module
cmake_minimum_required(VERSION 3.10)

project(CoFlame LANGUAGES Fortran VERSION 1.8)

find_package(MPI REQUIRED)

add_subdirectory(src)

install(PROGRAMS utils/tec2vtk DESTINATION bin)

