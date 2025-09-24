#include <pybind11/pybind11.h>
#include "../src/all.hpp"  // まとめヘッダを読み込む

namespace py = pybind11;

PYBIND11_MODULE(my_module, m) {
    m.doc() = "C++ functions exposed to Python via pybind11";

    // 例：a.cpp に int add(int, int) がある場合
    m.def("add", &add, "A function that adds two numbers");

    // 他の関数も同様に登録
    m.def("foo", &foo, "Description of foo");
    m.def("bar", &bar, "Description of bar");
}
