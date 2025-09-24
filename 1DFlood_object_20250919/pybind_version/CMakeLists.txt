cmake_minimum_required(VERSION 3.14)
project(mylib)

set(CMAKE_CXX_STANDARD 17)

set(pybind11_DIR "/home/tsuyoshi/my_env/lib/python3.12/site-packages/pybind11")

find_package(pybind11 REQUIRED)

pybind11_add_module(mylib
    src/ElemNod.cpp
    src/ReadWrite.cpp
    src/Solver.cpp
    python_bindings/bindings.cpp
)

pybind11_add_module(my_module ${SOURCES})