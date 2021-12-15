#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/contact_status.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_status, m) {
  py::enum_<ContactType>(m, "ContactType", py::arithmetic())
    .value("PointContact", ContactType::PointContact)
    .value("SurfaceContact", ContactType::SurfaceContact)
    .export_values();

  py::class_<ContactStatus>(m, "ContactStatus")
    .def(py::init<const int, const int>(),
          py::arg("max_point_contacts"), py::arg("contact_id")=0)
    .def("max_point_contacts", &ContactStatus::maxPointContacts)
    .def("is_contact_active", 
          static_cast<bool (ContactStatus::*)(const int) const>(&ContactStatus::isContactActive),
          py::arg("contact_index"))
    .def("is_contact_active", 
          static_cast<const std::vector<bool>& (ContactStatus::*)() const>(&ContactStatus::isContactActive))
    .def("activate_contact", &ContactStatus::activateContact,
          py::arg("contact_index"))
    .def("deactivate_contact", &ContactStatus::deactivateContact,
          py::arg("contact_index"))
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>& contact_indices)>(&ContactStatus::activateContacts))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>& contact_indices)>(&ContactStatus::deactivateContacts))
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)()>(&ContactStatus::activateContacts))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)()>(&ContactStatus::deactivateContacts))
    .def("set_contact_point", &ContactStatus::setContactPoint,
          py::arg("contact_index"), py::arg("contact_point"))
    .def("set_contact_points", &ContactStatus::setContactPoints,
          py::arg("contact_points"))
    .def("contact_point", &ContactStatus::contactPoint,
          py::arg("contact_index"))
    .def("contact_points", &ContactStatus::contactPoints)
    .def("set_contact_surface_rotation", &ContactStatus::setContactSurfaceRotation,
          py::arg("contact_index"), py::arg("contact_surface_rotation"))
    .def("set_contact_surfaces_rotations", &ContactStatus::setContactSurfacesRotations,
          py::arg("contact_surfaces_rotations"))
    .def("contact_surface_rotation", &ContactStatus::contactSurfaceRotation,
          py::arg("contact_index"))
    .def("contact_surfaces_rotations", &ContactStatus::contactSurfacesRotations)
    .def("set_contact_id", &ContactStatus::setContactId,
          py::arg("contact_id"))
    .def("contact_id", &ContactStatus::contactId)
    .def("__str__", [](const ContactStatus& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc