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
    .def(py::init<const int, const int>(),
          py::arg("max_point_contacts"), py::arg("impulse_mode_id")=0)
    .def("max_point_contacts", &ImpulseStatus::maxPointContacts)
    .def("is_impulse_active", 
          static_cast<bool (ImpulseStatus::*)(const int) const>(&ImpulseStatus::isImpulseActive),
          py::arg("contact_index"))
    .def("is_impulse_active", 
          static_cast<const std::vector<bool>& (ImpulseStatus::*)() const>(&ImpulseStatus::isImpulseActive))
    .def("activate_impulse", &ImpulseStatus::activateImpulse,
          py::arg("contact_index"))
    .def("deactivate_impulse", &ImpulseStatus::deactivateImpulse,
          py::arg("contact_index"))
    .def("activate_impulses", 
          static_cast<void (ImpulseStatus::*)(const std::vector<int>& impulse_indices)>(&ImpulseStatus::activateImpulses))
    .def("deactivate_impulses", 
          static_cast<void (ImpulseStatus::*)(const std::vector<int>& impulse_indices)>(&ImpulseStatus::deactivateImpulses))
    .def("activate_impulses", 
          static_cast<void (ImpulseStatus::*)()>(&ImpulseStatus::activateImpulses))
    .def("deactivate_impulses", 
          static_cast<void (ImpulseStatus::*)()>(&ImpulseStatus::deactivateImpulses))
    .def("set_contact_point", &ImpulseStatus::setContactPoint,
          py::arg("contact_index"), py::arg("contact_point"))
    .def("set_contact_points", &ImpulseStatus::setContactPoints,
          py::arg("contact_points"))
    .def("contact_point", &ImpulseStatus::contactPoint,
          py::arg("contact_index"))
    .def("contact_points", &ImpulseStatus::contactPoints)
    .def("set_contact_surface_rotation", &ImpulseStatus::setContactSurfaceRotation,
          py::arg("contact_index"), py::arg("contact_surface_rotation"))
    .def("set_contact_surfaces_rotations", &ImpulseStatus::setContactSurfacesRotations,
          py::arg("contact_surfaces_rotations"))
    .def("contact_surface_rotation", &ImpulseStatus::contactSurfaceRotation,
          py::arg("contact_index"))
    .def("contact_surfaces_rotations", &ImpulseStatus::contactSurfacesRotations)
    .def("set_impulse_mode_id", &ImpulseStatus::setImpulseId,
          py::arg("impulse_mode_id"))
    .def("impulse_mode_id", &ImpulseStatus::impulseId)
    .def("__str__", [](const ImpulseStatus& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc