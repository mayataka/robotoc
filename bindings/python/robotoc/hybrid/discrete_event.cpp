#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/discrete_event.hpp"

#include <iostream>


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(discrete_event, m) {
  py::enum_<DiscreteEventType>(m, "DiscreteEventType", py::arithmetic())
    .value("Impulse",  DiscreteEventType::Impulse)
    .value("Lift", DiscreteEventType::Lift)
    .value("NoneEvent", DiscreteEventType::None)
    .export_values();

  py::class_<DiscreteEvent>(m, "DiscreteEvent")
    .def(py::init<const std::vector<ContactType>&>(),
          py::arg("contact_types"))
    .def(py::init<const ContactStatus&, const ContactStatus&>(),
          py::arg("pre_contact_status"), py::arg("post_contact_status"))
    .def("clone", [](const DiscreteEvent& self) {
       auto other = self;
       return other;
     })
    .def("exist_discrete_event", &DiscreteEvent::existDiscreteEvent)
    .def("exist_impulse", &DiscreteEvent::existImpulse)
    .def("exist_lift", &DiscreteEvent::existLift)
    .def("impulse_status", &DiscreteEvent::impulseStatus)
    .def("pre_contact_status", &DiscreteEvent::preContactStatus)
    .def("post_contact_status", &DiscreteEvent::postContactStatus)
    .def("set_discrete_event", &DiscreteEvent::setDiscreteEvent,
          py::arg("pre_contact_status"), py::arg("post_contact_status"))
    .def("set_contact_placement", 
          static_cast<void (DiscreteEvent::*)(const int, const Eigen::Vector3d&)>(&DiscreteEvent::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"))
    .def("set_contact_placement", 
          static_cast<void (DiscreteEvent::*)(const int, const Eigen::Vector3d&, const Eigen::Matrix3d&)>(&DiscreteEvent::setContactPlacement),
          py::arg("contact_index"), py::arg("contact_position"), py::arg("contact_rotation"))
    .def("set_contact_placements", 
          static_cast<void (DiscreteEvent::*)(const std::vector<Eigen::Vector3d>&)>(&DiscreteEvent::setContactPlacements),
          py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (DiscreteEvent::*)(const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&DiscreteEvent::setContactPlacements),
          py::arg("contact_positions"), py::arg("contact_rotations"))
    .def("max_num_contacts", &DiscreteEvent::maxNumContacts)
    .def("event_type", &DiscreteEvent::eventType)
    .def("__str__", [](const DiscreteEvent& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc