#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/robot/impact_status.hpp"
#include "robotoc/utils/pybind11_macros.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(impact_status, m) {
  py::class_<ImpactStatus>(m, "ImpactStatus")
    .def(py::init<const std::vector<ContactType>&, const std::vector<std::string>&, const double, const int>(),
          py::arg("contact_types"), py::arg("contact_frame_names")=std::vector<std::string>({}), 
          py::arg("default_friction_coefficient")=0.7, py::arg("contact_mode_id")=0)
    .def(py::init<>())
    .def("max_num_contacts", &ImpactStatus::maxNumContacts)
    .def("is_impact_active", 
          static_cast<bool (ImpactStatus::*)(const int) const>(&ImpactStatus::isImpactActive),
          py::arg("contact_index"))
    .def("is_impact_active", 
          static_cast<const std::vector<bool>& (ImpactStatus::*)() const>(&ImpactStatus::isImpactActive))
    .def("activate_impact", &ImpactStatus::activateImpact,
          py::arg("contact_index"))
    .def("deactivate_impact", &ImpactStatus::deactivateImpact,
          py::arg("contact_index"))
    .def("activate_impacts", &ImpactStatus::activateImpacts,
          py::arg("contact_indices"))
    .def("deactivate_impacts", &ImpactStatus::deactivateImpacts,
          py::arg("contact_indices"))
    .def("set_contact_placement", 
          static_cast<void (ImpactStatus::*)(const int, const Eigen::Vector3d&)>(&ImpactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (ImpactStatus::*)(const int, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&ImpactStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placements", 
          static_cast<void (ImpactStatus::*)(const std::vector<Eigen::Vector3d>&)>(&ImpactStatus::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ImpactStatus::*)(const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&ImpactStatus::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("contact_placement", &ImpactStatus::contactPlacement,
          py::arg("contact_index"))
    .def("contact_placements", &ImpactStatus::contactPlacements)
    .def("contact_position", &ImpactStatus::contactPosition,
          py::arg("contact_index"))
    .def("contact_positions", &ImpactStatus::contactPositions)
    .def("contact_rotation", &ImpactStatus::contactRotation,
          py::arg("contact_index"))
    .def("contact_rotations", &ImpactStatus::contactRotations)
    .def("set_friction_coefficient", 
          static_cast<void (ImpactStatus::*)(const int, const double)>(&ImpactStatus::setFrictionCoefficient),
          py::arg("contact_index"), py::arg("friction_coefficient"))
    .def("set_friction_coefficient", 
          static_cast<void (ImpactStatus::*)(const std::string&, const double)>(&ImpactStatus::setFrictionCoefficient),
          py::arg("contact_frame_name"), py::arg("friction_coefficient"))
    .def("set_friction_coefficients", 
          static_cast<void (ImpactStatus::*)(const std::vector<double>&)>(&ImpactStatus::setFrictionCoefficients),
          py::arg("friction_coefficients"))
    .def("set_friction_coefficients", 
          static_cast<void (ImpactStatus::*)(const std::unordered_map<std::string, double>&)>(&ImpactStatus::setFrictionCoefficients),
          py::arg("friction_coefficients"))
    .def("friction_coefficient", 
          static_cast<double (ImpactStatus::*)(const int) const>(&ImpactStatus::frictionCoefficient),
          py::arg("contact_index"))
    .def("friction_coefficient", 
          static_cast<double (ImpactStatus::*)(const std::string&) const>(&ImpactStatus::frictionCoefficient),
          py::arg("contact_frame_name"))
    .def("friction_coefficients", &ImpactStatus::frictionCoefficients)
    .def("set_impact_mode_id", &ImpactStatus::setImpactModeId,
          py::arg("impact_mode_id"))
    .def("impact_mode_id", &ImpactStatus::impactModeId)
     DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ImpactStatus)
     DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(ImpactStatus);
}

} // namespace python
} // namespace robotoc