#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/ocp/split_solution.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_solution, m) {
  py::class_<SplitSolution>(m, "SplitSolution")
    .def(py::init<const Robot&>())
    .def_readwrite("q", &SplitSolution::q)
    .def_readwrite("v", &SplitSolution::v)
    .def_readwrite("a", &SplitSolution::a)
    .def_readwrite("u", &SplitSolution::u)
    .def_readwrite("f", &SplitSolution::f)
    .def_readwrite("lmd", &SplitSolution::lmd)
    .def_readwrite("gmm", &SplitSolution::gmm)
    .def_readwrite("beta", &SplitSolution::beta)
    .def_readwrite("mu", &SplitSolution::mu)
    .def("f_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::f_stack))
    .def("mu_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::mu_stack))
    .def("xi_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::xi_stack))
    .def("__str__", [](const SplitSolution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc