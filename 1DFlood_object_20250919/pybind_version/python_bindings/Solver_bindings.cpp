#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "../src/all.hpp"
namespace py = pybind11;

void bind_Time_solver(py::module_ &m) {
    py::class_<Time_solver, std::shared_ptr<Time_solver>>(m, "Time_solver")
        .def("update_depth", &Time_solver::update_depth) ;
    py::class_<Euler, Time_solver, std::shared_ptr<Euler>>(m,"Euler")
        .def(py::init<>()) ;
    py::class_<ButcherTable>(m, "ButcherTable") 
        .def_readonly("stages", &ButcherTable::stages)
        .def_readonly("stage_weights", &ButcherTable::stage_weights)
        .def_readonly("final_weights", &ButcherTable::final_weights) ;
    py::class_<Runge_Kutta, Time_solver, std::shared_ptr<Runge_Kutta>>(m,"Runge_Kutta")
        .def(py::init<const ButcherTable&>()) ;
        m.def("RK2", &RKTables::RK2, py::return_value_policy::reference);
        m.def("RK3", &RKTables::RK3, py::return_value_policy::reference);
        m.def("RK4", &RKTables::RK4, py::return_value_policy::reference);
        m.def("RK6", &RKTables::RK6, py::return_value_policy::reference);        
    m.def("compute_all", &compute_all) ;
}
