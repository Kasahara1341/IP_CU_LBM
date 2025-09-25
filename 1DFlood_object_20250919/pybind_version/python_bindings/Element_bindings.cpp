#include <pybind11/pybind11.h>
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

void bind_Node(py::module_ &m) {
    py::class_<Node, std::shared_ptr<Node>>(m, "Node")
        .def(py::init<>())
        // .def_readwrite("position", &Element::position)
        .def("solve_momentum_equation", &Node::solve_momentum_equation)
        .def("set_flux", &Node::set_flow_q)
        .def("set_up_element", &Node::set_up_element)
        .def("set_dn_element", &Node::set_dn_element)

        .def("get_flux",       &Node::get_flux)
        .def("get_up_element", &Node::get_up_element)
        .def("get_dn_element", &Node::get_dn_element) ;
}

void bind_Time_solver(py::module_ &m) {
    py::class_<Time_solver, std::shared_ptr<Time_solver>>(m, "Time_solver")
        .def("update_depth", &Time_solver::update_depth) ;
    py::class_<Euler, Time_solver, std::shared_ptr<Euler>>(m,"Euler")
        .def(py::init<>()) ;
    py::class_<ButcherTable>(m, "ButcherTable")
        .def_readonly("stages", &ButcherTable::stages)
        .def_readonly("stage_weights", &ButcherTable::stage_weights)
        .def_readonly("final_weights", &ButcherTable::final_weights) ;
        m.def("RK2", &RKTables::RK2, py::return_value_policy::reference);
        m.def("RK3", &RKTables::RK3, py::return_value_policy::reference);
        m.def("RK4", &RKTables::RK4, py::return_value_policy::reference);
        m.def("RK6", &RKTables::RK6, py::return_value_policy::reference);        
    py::class_<Runge_Kutta, Time_solver, std::shared_ptr<Runge_Kutta>>(m,"Runge_Kutta")
        .def(py::init<const ButcherTable&>()) ;
}