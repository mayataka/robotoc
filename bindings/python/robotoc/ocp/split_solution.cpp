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
    .def(py::init<>())
    .def("clone", [](const SplitSolution& self) {
       auto other = self;
       return other;
     })
    .def("set_contact_status", 
          static_cast<void (SplitSolution::*)(const ContactStatus&)>(&SplitSolution::setContactStatus),
          py::arg("contact_status"))
    .def("set_contact_status", 
          static_cast<void (SplitSolution::*)(const SplitSolution&)>(&SplitSolution::setContactStatus),
          py::arg("other"))
    .def("set_impulse_status", 
          static_cast<void (SplitSolution::*)(const ImpulseStatus&)>(&SplitSolution::setImpulseStatus),
          py::arg("impulse_status"))
    .def("set_impulse_status", 
          static_cast<void (SplitSolution::*)(const SplitSolution&)>(&SplitSolution::setImpulseStatus),
          py::arg("other"))
    .def("set_impulse_status", 
          static_cast<void (SplitSolution::*)()>(&SplitSolution::setImpulseStatus))
    .def("set_f_stack", &SplitSolution::set_f_stack)
    .def("set_f_vector", &SplitSolution::set_f_vector)
    .def("set_mu_stack", &SplitSolution::set_mu_stack)
    .def("set_mu_vector", &SplitSolution::set_mu_vector)
    .def("dimf", &SplitSolution::dimf)
    .def("dimi", &SplitSolution::dimi)
    .def("has_active_contacts", &SplitSolution::hasActiveContacts)
    .def("has_active_impulse", &SplitSolution::hasActiveImpulse)
    .def("is_contact_active", 
          static_cast<bool (SplitSolution::*)(const int) const>(&SplitSolution::isContactActive),
          py::arg("contact_index"))
    .def("is_contact_active", 
          static_cast<std::vector<bool> (SplitSolution::*)() const>(&SplitSolution::isContactActive))
    .def_readwrite("q", &SplitSolution::q)
    .def_readwrite("v", &SplitSolution::v)
    .def_readwrite("a", &SplitSolution::a)
    .def_readwrite("u", &SplitSolution::u)
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
    .def("__str__", [](const SplitSolution& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc