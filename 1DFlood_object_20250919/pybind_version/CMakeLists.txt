cmake_minimum_required(VERSION 3.14)
project(my_module)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# option
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# pybind11 を FetchContent で取得
include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.13.1  # 適宜最新版に変更
)

FetchContent_MakeAvailable(pybind11)

# Python モジュールとしてビルド
pybind11_add_module(my_module 
    python_bindings/bindings.cpp
    python_bindings/Element_bindings.cpp
    python_bindings/Node_bindings.cpp
    python_bindings/Solver_bindings.cpp
    src/Element.cpp
    src/Node.cpp
    src/Solver.cpp
)
