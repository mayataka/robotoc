#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/impulse_status.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impulse_status, m) {
  py::class_<ImpulseStatus>(m, "ImpulseStatus")
    .def(py::init<const int>())
    .def("max_point_contacts", &ImpulseStatus::maxPointContacts)
    .def("is_impulse_active", 
          static_cast<bool (ImpulseStatus::*)(const int) const>(&ImpulseStatus::isImpulseActive))
    .def("is_impulse_active", 
          static_cast<const std::vector<bool>& (ImpulseStatus::*)() const>(&ImpulseStatus::isImpulseActive))
    .def("activate_impulse", &ImpulseStatus::activateImpulse)
    .def("deactivate_impulse", &ImpulseStatus::deactivateImpulse)
    .def("activate_impulses", 
          static_cast<void (ImpulseStatus::*)(const std::vector<int>& impulse_indices)>(&ImpulseStatus::activateImpulses))
    .def("deactivate_impulses", 
          static_cast<void (ImpulseStatus::*)(const std::vector<int>& impulse_indices)>(&ImpulseStatus::deactivateImpulses))
    .def("activate_impulses", 
          static_cast<void (ImpulseStatus::*)()>(&ImpulseStatus::activateImpulses))
    .def("deactivate_impulses", 
          static_cast<void (ImpulseStatus::*)()>(&ImpulseStatus::deactivateImpulses))
    .def("set_contact_point", &ImpulseStatus::setContactPoint)
    .def("set_contact_points", &ImpulseStatus::setContactPoints)
    .def("contact_point", &ImpulseStatus::contactPoint)
    .def("contact_points", &ImpulseStatus::contactPoints)
    .def("__str__", [](const ImpulseStatus& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc