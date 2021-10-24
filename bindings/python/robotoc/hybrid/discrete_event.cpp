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
    .def(py::init<const int>())
    .def(py::init<const ContactStatus&, const ContactStatus&>())
    .def("exist_discrete_event", &DiscreteEvent::existDiscreteEvent)
    .def("exist_impulse", &DiscreteEvent::existImpulse)
    .def("exist_lift", &DiscreteEvent::existLift)
    .def("impulse_status", &DiscreteEvent::impulseStatus)
    .def("pre_contact_status", &DiscreteEvent::preContactStatus)
    .def("post_contact_status", &DiscreteEvent::postContactStatus)
    .def("set_discrete_event", &DiscreteEvent::setDiscreteEvent)
    .def("set_contact_point", &DiscreteEvent::setContactPoint,
          py::arg("contact_index"), py::arg("contact_point"))
    .def("set_contact_points", &DiscreteEvent::setContactPoints,
          py::arg("contact_points"))
    .def("max_point_contacts", &DiscreteEvent::maxPointContacts)
    .def("event_type", &DiscreteEvent::eventType)
    .def("__str__", [](const DiscreteEvent& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc