#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/core/split_solution.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(split_solution, m) {
  py::class_<SplitSolution>(m, "SplitSolution")
    .def(py::init<const Robot&>())
    .def(py::init<>())
    .def("set_contact_status", 
          static_cast<void (SplitSolution::*)(const ContactStatus&)>(&SplitSolution::setContactStatus),
          py::arg("contact_status"))
    .def("set_contact_status", 
          static_cast<void (SplitSolution::*)(const ImpactStatus&)>(&SplitSolution::setContactStatus),
          py::arg("contact_status"))
    .def("set_contact_status", 
          static_cast<void (SplitSolution::*)(const SplitSolution&)>(&SplitSolution::setContactStatus),
          py::arg("other"))
    .def("set_switching_constraint_dimension", &SplitSolution::setSwitchingConstraintDimension,
          py::arg("dims"))
    .def("set_f_stack", &SplitSolution::set_f_stack)
    .def("set_f_vector", &SplitSolution::set_f_vector)
    .def("set_mu_stack", &SplitSolution::set_mu_stack)
    .def("set_mu_vector", &SplitSolution::set_mu_vector)
    .def("dimf", &SplitSolution::dimf)
    .def("dims", &SplitSolution::dims)
    .def("is_contact_active", 
          static_cast<bool (SplitSolution::*)(const int) const>(&SplitSolution::isContactActive),
          py::arg("contact_index"))
    .def("is_contact_active", 
          static_cast<std::vector<bool> (SplitSolution::*)() const>(&SplitSolution::isContactActive))
    .def_readwrite("q", &SplitSolution::q)
    .def_readwrite("v", &SplitSolution::v)
    .def_readwrite("u", &SplitSolution::u)
    .def_readwrite("a", &SplitSolution::a)
    .def_readwrite("dv", &SplitSolution::dv)
    .def_readwrite("f", &SplitSolution::f)
    .def_readwrite("lmd", &SplitSolution::lmd)
    .def_readwrite("gmm", &SplitSolution::gmm)
    .def_readwrite("beta", &SplitSolution::beta)
    .def_readwrite("mu", &SplitSolution::mu)
    .def_readwrite("nu_passive", &SplitSolution::nu_passive)
    .def_property("f_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::f_stack),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitSolution::*)()>(&SplitSolution::f_stack))
    .def_property("mu_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::mu_stack),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitSolution::*)()>(&SplitSolution::mu_stack))
    .def_property("xi_stack", static_cast<const Eigen::VectorBlock<const Eigen::VectorXd> (SplitSolution::*)() const>(&SplitSolution::xi_stack),
                             static_cast<Eigen::VectorBlock<Eigen::VectorXd> (SplitSolution::*)()>(&SplitSolution::xi_stack))
    .def("dimf", &SplitSolution::dimf)
    .def("dims", &SplitSolution::dims)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(SplitSolution)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(SplitSolution);
}

} // namespace python
} // namespace robotoc