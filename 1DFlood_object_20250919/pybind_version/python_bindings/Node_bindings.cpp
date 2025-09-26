#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "../src/all.hpp"
namespace py = pybind11;

void bind_Node(py::module_ &m) {
    py::class_<Node, std::shared_ptr<Node>>(m, "Node")
        .def(py::init<>())
        // .def_readwrite("position", &Element::position)
        .def("solve_momentum_equation", &Node::solve_momentum_equation)
        .def("set_flux", &Node::set_flux)
        .def("set_up_element", &Node::set_up_element)
        .def("set_dn_element", &Node::set_dn_element)

        .def("get_flux",       &Node::get_flux)
        .def("get_up_element", &Node::get_up_element)
        .def("get_dn_element", &Node::get_dn_element) ;
}