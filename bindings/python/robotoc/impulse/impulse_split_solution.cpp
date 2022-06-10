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
    .def(py::init<>())
    .def("set_impulse_status", 
          static_cast<void (ImpulseSplitSolution::*)(const ImpulseStatus&)>(&ImpulseSplitSolution::setImpulseStatus),
          py::arg("impulse_status"))
    .def("set_impulse_status", 
          static_cast<void (ImpulseSplitSolution::*)(const ImpulseSplitSolution&)>(&ImpulseSplitSolution::setImpulseStatus),
          py::arg("other"))
    .def("set_f_stack", &ImpulseSplitSolution::set_f_stack)
    .def("set_f_vector", &ImpulseSplitSolution::set_f_vector)
    .def("set_mu_stack", &ImpulseSplitSolution::set_mu_stack)
    .def("set_mu_vector", &ImpulseSplitSolution::set_mu_vector)
    .def("dimi", &ImpulseSplitSolution::dimi)
    .def("is_impulse_active", 
          static_cast<bool (ImpulseSplitSolution::*)(const int) const>(&ImpulseSplitSolution::isImpulseActive),
          py::arg("contact_index"))
    .def("is_impulse_active", 
          static_cast<std::vector<bool> (ImpulseSplitSolution::*)() const>(&ImpulseSplitSolution::isImpulseActive))
    .def_readwrite("q", &ImpulseSplitSolution::q)
    .def_readwrite("v", &ImpulseSplitSolution::v)
    .def_readwrite("dv", &ImpulseSplitSolution::dv)
    .def_readwrite("f", &ImpulseSplitSolution::f)
    .def_readwrite("lmd", &ImpulseSplitSolution::lmd)
    .def_readwrite("gmm", &ImpulseSplitSolution::gmm)
    .def_readwrite("beta", &ImpulseSplitSolution::beta)
    .def_readwrite("mu", &ImpulseSplitSolution::mu)
    .def_property("f_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitSolution::*)() const>(&ImpulseSplitSolution::f_stack),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitSolution::*)()>(&ImpulseSplitSolution::f_stack))
    .def_property("mu_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (ImpulseSplitSolution::*)() const>(&ImpulseSplitSolution::mu_stack),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (ImpulseSplitSolution::*)()>(&ImpulseSplitSolution::mu_stack))
    .def("__str__", [](const ImpulseSplitSolution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc