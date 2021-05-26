#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/ocp/split_solution.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_solution, m) {
  py::class_<SplitSolution>(m, "SplitSolution")
    .def(py::init<const Robot&>())
    .def_readwrite("q", &SplitSolution::q)
    .def_readwrite("v", &SplitSolution::v)
    .def_readwrite("a", &SplitSolution::a)
    .def_readwrite("u", &SplitSolution::u)
    .def_readwrite("f", &SplitSolution::f);
}

} // namespace python
} // namespace idocp