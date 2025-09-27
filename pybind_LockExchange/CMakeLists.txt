cmake_minimum_required(VERSION 3.14)
project(cuda_pybind_example LANGUAGES CXX CUDA)

find_package(pybind11 REQUIRED)

add_library(my_module SHARED wrapper.cpp kernel.cu)
set_target_properties(my_module PROPERTIES PREFIX "")

target_compile_features(my_module PUBLIC cxx_std_17)
target_compile_options(my_module PRIVATE 
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr -arch=sm_89>)
    
find_package(CUDAToolkit REQUIRED)

target_link_libraries(my_module PRIVATE 
    pybind11::module 
    CUDA::cudart)
