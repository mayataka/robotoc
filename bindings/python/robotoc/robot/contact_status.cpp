#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/contact_status.hpp"
#include "robotoc/utils/pybind11_macros.hpp"

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
    .def(py::init<const std::vector<ContactType>&, const std::vector<std::string>&, const double>(),
          py::arg("contact_types"), py::arg("contact_frame_names")=std::vector<std::string>({}), 
          py::arg("default_friction_coefficient")=0.7)
    .def(py::init<>())
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
    .def("set_friction_coefficient", 
          static_cast<void (ContactStatus::*)(const int, const double)>(&ContactStatus::setFrictionCoefficient),
          py::arg("contact_index"), py::arg("friction_coefficient"))
    .def("set_friction_coefficient", 
          static_cast<void (ContactStatus::*)(const std::string&, const double)>(&ContactStatus::setFrictionCoefficient),
          py::arg("contact_frame_name"), py::arg("friction_coefficient"))
    .def("set_friction_coefficients", 
          static_cast<void (ContactStatus::*)(const std::vector<double>&)>(&ContactStatus::setFrictionCoefficients),
          py::arg("friction_coefficients"))
    .def("set_friction_coefficients", 
          static_cast<void (ContactStatus::*)(const std::unordered_map<std::string, double>&)>(&ContactStatus::setFrictionCoefficients),
          py::arg("friction_coefficients"))
    .def("friction_coefficient", 
          static_cast<double (ContactStatus::*)(const int) const>(&ContactStatus::frictionCoefficient),
          py::arg("contact_index"))
    .def("friction_coefficient", 
          static_cast<double (ContactStatus::*)(const std::string&) const>(&ContactStatus::frictionCoefficient),
          py::arg("contact_frame_name"))
    .def("friction_coefficients", &ContactStatus::frictionCoefficients)
    .def("find_contact_index", &ContactStatus::findContactIndex,
          py::arg("contact_frame_name"))
     DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ContactStatus)
     DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(ContactStatus);
}

} // namespace python
} // namespace robotoc