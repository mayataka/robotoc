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
    .def(py::init<const std::vector<ContactType>&, const int>(),
          py::arg("contact_types"), py::arg("contact_mode_id")=0)
    .def("max_num_contacts", &ContactStatus::maxNumContacts)
    .def("is_contact_active", 
          static_cast<bool (ContactStatus::*)(const int) const>(&ContactStatus::isContactActive),
          py::arg("contact_index"))
    .def("is_contact_active", 
          static_cast<const std::vector<bool>& (ContactStatus::*)() const>(&ContactStatus::isContactActive))
    .def("activate_contact", &ContactStatus::activateContact,
          py::arg("contact_index"))
    .def("deactivate_contact", &ContactStatus::deactivateContact,
          py::arg("contact_index"))
    .def("activate_contacts", &ContactStatus::activateContacts,
          py::arg("contact_indices"))
    .def("deactivate_contacts", &ContactStatus::deactivateContacts,
          py::arg("contact_indices"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const int, const Eigen::Vector3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const int, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::vector<Eigen::Vector3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("contact_placement", &ContactStatus::contactPlacement,
          py::arg("contact_index"))
    .def("contact_placements", &ContactStatus::contactPlacements)
    .def("contact_position", &ContactStatus::contactPosition,
          py::arg("contact_index"))
    .def("contact_positions", &ContactStatus::contactPositions)
    .def("contact_rotation", &ContactStatus::contactRotation,
          py::arg("contact_index"))
    .def("contact_rotations", &ContactStatus::contactRotations)
    .def("set_contact_mode_id", &ContactStatus::setContactModeId,
          py::arg("contact_mode_id"))
    .def("contact_mode_id", &ContactStatus::contactModeId)
    .def("__str__", [](const ContactStatus& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc