#include <pybind11/pybind11.h>
#include "../src/all.hpp"

namespace py = pybind11;

// ""内のものはpythonファイルで呼ぶときの名前でここで決められる
void bind_Element(py::module_ &m) ;
void bind_Node(py::module_ &m) ;
void bind_Time_solver(py::module_ &m) ;

//  (pythonでのmodule名  ここでの呼び名)
PYBIND11_MODULE(my_module, m) {
    bind_Element(m);
    bind_Node(m);
    bind_Time_solver(m);
}
