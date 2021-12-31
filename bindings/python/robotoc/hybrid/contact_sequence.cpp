#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/hybrid/contact_sequence.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_sequence, m) {
  py::class_<ContactSequence, std::shared_ptr<ContactSequence>>(m, "ContactSequence")
    .def(py::init<const Robot&, const int>(),
         py::arg("robot"), py::arg("max_num_each_events")=0)
    .def("init_contact_sequence", &ContactSequence::initContactSequence,
           py::arg("contact_status"))
    .def("push_back", static_cast<void (ContactSequence::*)(const DiscreteEvent&, const double, const bool)>(&ContactSequence::push_back),
          py::arg("discrete_event"), py::arg("event_time"), py::arg("sto")=false)
    .def("push_back", static_cast<void (ContactSequence::*)(const ContactStatus&, const double, const bool)>(&ContactSequence::push_back),
          py::arg("contact_status"), py::arg("switching_time"), py::arg("sto")=false)
    .def("pop_back", &ContactSequence::pop_back)
    .def("pop_front", &ContactSequence::pop_front)
    .def("set_impulse_time", &ContactSequence::setImpulseTime, 
          py::arg("impulse_index"), py::arg("impulse_time"))
    .def("set_lift_time", &ContactSequence::setLiftTime, 
          py::arg("lift_index"), py::arg("lift_time"))
    .def("is_STO_enabled_impulse", &ContactSequence::isSTOEnabledImpulse, 
          py::arg("impulse_index"))
    .def("is_STO_enabled_lift", &ContactSequence::isSTOEnabledLift, 
          py::arg("lift_index"))
    .def("is_event_time_consistent", &ContactSequence::isEventTimeConsistent)
    .def("set_contact_placements", 
          static_cast<void (ContactSequence::*)(const int, const std::vector<Eigen::Vector3d>&)>(&ContactSequence::setContactPlacements),
          py::arg("contact_phase"), py::arg("contact_positions"))
    .def("set_contact_placements", 
          static_cast<void (ContactSequence::*)(const int, const std::vector<Eigen::Vector3d>&, const std::vector<Eigen::Matrix3d>&)>(&ContactSequence::setContactPlacements),
          py::arg("contact_phase"), py::arg("contact_positions"), 
          py::arg("contact_rotations"))
    .def("num_impulse_events", &ContactSequence::numImpulseEvents)
    .def("num_lift_events", &ContactSequence::numLiftEvents)
    .def("num_discrete_events", &ContactSequence::numDiscreteEvents)
    .def("num_contact_phases", &ContactSequence::numContactPhases)
    .def("contact_status", &ContactSequence::contactStatus,
          py::arg("contact_phase"))
    .def("impulse_status", &ContactSequence::impulseStatus,
          py::arg("impulse_index"))
    .def("impulse_time", &ContactSequence::impulseTime,
          py::arg("impulse_index"))
    .def("lift_time", &ContactSequence::liftTime,
          py::arg("lift_index"))
    .def("event_type", &ContactSequence::eventType,
          py::arg("event_index"))
    .def("event_times", &ContactSequence::eventTimes)
    .def("max_num_each_events", &ContactSequence::maxNumEachEvents)
    .def("max_num_events", &ContactSequence::maxNumEvents)
    .def("__str__", [](const ContactSequence& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc