#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/impulse/impulse_split_solution.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_split_solution, m) {
  py::class_<ImpulseSplitSolution>(m, "ImpulseSplitSolution")
    .def(py::init<const Robot&>())
    .def_readwrite("q", &ImpulseSplitSolution::q)
    .def_readwrite("v", &ImpulseSplitSolution::v)
    .def_readwrite("dv", &ImpulseSplitSolution::dv)
    .def_readwrite("f", &ImpulseSplitSolution::f)
    .def_readwrite("lmd", &ImpulseSplitSolution::lmd)
    .def_readwrite("gmm", &ImpulseSplitSolution::gmm)
    .def_readwrite("beta", &ImpulseSplitSolution::beta)
    .def_readwrite("mu", &ImpulseSplitSolution::mu)
    .def("f_stack", 
         static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitSolution::*)() const>(&ImpulseSplitSolution::f_stack))
    .def("mu_stack", 
         static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitSolution::*)() const>(&ImpulseSplitSolution::mu_stack))
    .def("__str__", [](const ImpulseSplitSolution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc