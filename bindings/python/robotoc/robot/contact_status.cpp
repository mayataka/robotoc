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
    .def(py::init<const std::vector<ContactType>&, const std::vector<std::string>&, const int>(),
          py::arg("contact_types"), py::arg("contact_frame_names")=std::vector<std::string>({}), 
          py::arg("contact_mode_id")=0)
    .def("clone", [](const ContactStatus& self) {
       auto other = self;
       return other;
     })
    .def("max_num_contacts", &ContactStatus::maxNumContacts)
    .def("is_contact_active", 
          static_cast<bool (ContactStatus::*)(const int) const>(&ContactStatus::isContactActive),
          py::arg("contact_index"))
    .def("is_contact_active", 
          static_cast<bool (ContactStatus::*)(const std::string&) const>(&ContactStatus::isContactActive),
          py::arg("contact_frame_name"))
    .def("is_contact_active", 
          static_cast<const std::vector<bool>& (ContactStatus::*)() const>(&ContactStatus::isContactActive))
    .def("activate_contact", 
          static_cast<void (ContactStatus::*)(const int)>(&ContactStatus::activateContact),
          py::arg("contact_index"))
    .def("activate_contact", 
          static_cast<void (ContactStatus::*)(const std::string&)>(&ContactStatus::activateContact),
          py::arg("contact_frame_name"))
    .def("deactivate_contact", 
          static_cast<void (ContactStatus::*)(const int)>(&ContactStatus::deactivateContact),
          py::arg("contact_index"))
    .def("deactivate_contact", 
          static_cast<void (ContactStatus::*)(const std::string&)>(&ContactStatus::deactivateContact),
          py::arg("contact_frame_name"))
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>&)>(&ContactStatus::activateContacts),
          py::arg("contact_indices"))
    .def("activate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<std::string>&)>(&ContactStatus::activateContacts),
          py::arg("contact_frame_names"))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<int>&)>(&ContactStatus::deactivateContacts),
          py::arg("contact_indices"))
    .def("deactivate_contacts", 
          static_cast<void (ContactStatus::*)(const std::vector<std::string>&)>(&ContactStatus::deactivateContacts),
          py::arg("contact_frame_names"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const int, const Eigen::Vector3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const std::string&, const Eigen::Vector3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_frame_name"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const int, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const std::string&, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_frame_name"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const int, const SE3&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_placement"))
    .def("set_contact_placement", 
          static_cast<void (ContactStatus::*)(const std::string&, const SE3&)>(&ContactStatus::setContactPlacement),
          py::arg("contact_frame_name"), py::arg("contact_placement"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::vector<Eigen::Vector3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::unordered_map<std::string, Eigen::Vector3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const std::unordered_map<std::string, Eigen::Vector3d>&, 
                                              const std::unordered_map<std::string, Eigen::Matrix3d>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const aligned_vector<SE3>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_placements"))
    .def("set_contact_placements", 
          static_cast<void (ContactStatus::*)(const aligned_unordered_map<std::string, SE3>&)>(&ContactStatus::setContactPlacements),
          py::arg("contact_placements"))
    .def("contact_placement", 
          static_cast<const SE3& (ContactStatus::*)(const int) const>(&ContactStatus::contactPlacement),
          py::arg("contact_index"))
    .def("contact_placement", 
          static_cast<const SE3& (ContactStatus::*)(const std::string&) const>(&ContactStatus::contactPlacement),
          py::arg("contact_frame_name"))
    .def("contact_placements", &ContactStatus::contactPlacements)
    .def("contact_position", 
          static_cast<const Eigen::Vector3d& (ContactStatus::*)(const int) const>(&ContactStatus::contactPosition),
          py::arg("contact_index"))
    .def("contact_position", 
          static_cast<const Eigen::Vector3d& (ContactStatus::*)(const std::string&) const>(&ContactStatus::contactPosition),
          py::arg("contact_frame_name"))
    .def("contact_positions", &ContactStatus::contactPositions)
    .def("contact_rotation", 
          static_cast<const Eigen::Matrix3d& (ContactStatus::*)(const int) const>(&ContactStatus::contactRotation),
          py::arg("contact_index"))
    .def("contact_rotation", 
          static_cast<const Eigen::Matrix3d& (ContactStatus::*)(const std::string&) const>(&ContactStatus::contactRotation),
          py::arg("contact_frame_name"))
    .def("contact_rotations", &ContactStatus::contactRotations)
    .def("find_contact_index", &ContactStatus::findContactIndex,
          py::arg("contact_frame_name"))
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