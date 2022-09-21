#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/solution.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(solution, m) {
  py::class_<Solution>(m, "Solution")
    .def(py::init<const Robot&, const int, const int>(),
         py::arg("robot"), py::arg("N"), py::arg("reserved_num_discrete_eventsbve"))
    .def(py::init<>())
    .def("__getitem__", [](const Solution& self, const int i) {
        return self.data[i];
      })
    .def("__setitem__", [](Solution& self, const int i, const SplitSolution& s) {
        self.data[i] = s;
      })
    .def_readwrite("aux", &Solution::aux)
    .def_readwrite("lift", &Solution::lift)
    .def_readwrite("impulse", &Solution::impulse)
    .def("__str__", [](const Solution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc