#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "../src/all.hpp"
namespace py = pybind11;

void bind_Element(py::module_ &m) {
    py::class_<Element, std::shared_ptr<Element>>(m, "Element")
        .def(py::init<double, double, double, double, double>())
        // .def_readwrite("position", &Element::position)
        .def("set_depth",           &Element::set_depth)
        .def("set_up_node",         &Element::set_up_node)
        .def("set_dn_node",         &Element::set_dn_node)
        .def("set_time_solver",     &Element::set_time_solver)
        .def("solve_mass_equation", &Element::solve_mass_equation)
        .def("calc_increment",      &Element::calc_increment)
        .def("get_elevation",       &Element::get_elevation)
        .def("get_depth",           &Element::get_depth)
        .def("get_length",          &Element::get_length)
        .def("get_width",           &Element::get_width)
        .def("get_position",        &Element::get_position) ;
        
}
