#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/impulse/impulse_split_solution.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_split_solution, m) {
  py::class_<ImpulseSplitSolution>(m, "ImpulseSplitSolution")
    .def(py::init<const Robot&>())
    .def_readwrite("q", &ImpulseSplitSolution::q)
    .def_readwrite("v", &ImpulseSplitSolution::v)
    .def_readwrite("dv", &ImpulseSplitSolution::dv)
    .def_readwrite("f", &ImpulseSplitSolution::f)
    .def("__str__", [](const ImpulseSplitSolution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace idocp