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
    .def(py::init<const std::vector<ContactType>&, const std::vector<std::string>&, const int>(),
          py::arg("contact_types"), py::arg("contact_frame_names")=std::vector<std::string>({}), 
          py::arg("contact_mode_id")=0)
    .def("clone", [](const ImpulseStatus& self) {
       auto other = self;
       return other;
     })
    .def("max_num_contacts", &ImpulseStatus::maxNumContacts)
    .def("is_impulse_active", 
          static_cast<bool (ImpulseStatus::*)(const int) const>(&ImpulseStatus::isImpulseActive),
          py::arg("contact_index"))
    .def("is_impulse_active", 
          static_cast<const std::vector<bool>& (ImpulseStatus::*)() const>(&ImpulseStatus::isImpulseActive))
    .def("activate_impulse", &ImpulseStatus::activateImpulse,
          py::arg("contact_index"))
    .def("deactivate_impulse", &ImpulseStatus::deactivateImpulse,
          py::arg("contact_index"))
    .def("activate_impulses", &ImpulseStatus::activateImpulses,
          py::arg("contact_indices"))
    .def("deactivate_impulses", &ImpulseStatus::deactivateImpulses,
          py::arg("contact_indices"))
    .def("set_contact_placement", 
          static_cast<void (ImpulseStatus::*)(const int, const Eigen::Vector3d&)>(&ImpulseStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (ImpulseStatus::*)(const int, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&ImpulseStatus::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placements", 
          static_cast<void (ImpulseStatus::*)(const std::vector<Eigen::Vector3d>&)>(&ImpulseStatus::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ImpulseStatus::*)(const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&ImpulseStatus::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("contact_placement", &ImpulseStatus::contactPlacement,
          py::arg("contact_index"))
    .def("contact_placements", &ImpulseStatus::contactPlacements)
    .def("contact_position", &ImpulseStatus::contactPosition,
          py::arg("contact_index"))
    .def("contact_positions", &ImpulseStatus::contactPositions)
    .def("contact_rotation", &ImpulseStatus::contactRotation,
          py::arg("contact_index"))
    .def("contact_rotations", &ImpulseStatus::contactRotations)
    .def("set_impulse_mode_id", &ImpulseStatus::setImpulseModeId,
          py::arg("impulse_mode_id"))
    .def("impulse_mode_id", &ImpulseStatus::impulseModeId)
    .def("__str__", [](const ImpulseStatus& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc