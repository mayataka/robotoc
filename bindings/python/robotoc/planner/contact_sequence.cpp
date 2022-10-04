#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(contact_sequence, m) {
  py::class_<ContactSequence, std::shared_ptr<ContactSequence>>(m, "ContactSequence")
    .def(py::init<const Robot&, const int>(),
         py::arg("robot"), py::arg("reserved_num_discrete_events")=0)
    .def("init", &ContactSequence::init,
           py::arg("contact_status"))
    .def("push_back", static_cast<void (ContactSequence::*)(const DiscreteEvent&, const double, const bool)>(&ContactSequence::push_back),
          py::arg("discrete_event"), py::arg("event_time"), py::arg("sto")=false)
    .def("push_back", static_cast<void (ContactSequence::*)(const ContactStatus&, const double, const bool)>(&ContactSequence::push_back),
          py::arg("contact_status"), py::arg("switching_time"), py::arg("sto")=false)
    .def("pop_back", &ContactSequence::pop_back)
    .def("pop_front", &ContactSequence::pop_front)
    .def("set_impact_time", &ContactSequence::setImpactTime, 
          py::arg("impact_index"), py::arg("impact_time"))
    .def("set_lift_time", &ContactSequence::setLiftTime, 
          py::arg("lift_index"), py::arg("lift_time"))
    .def("is_STO_enabled_impact", &ContactSequence::isSTOEnabledImpact, 
          py::arg("impact_index"))
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
    .def("set_friction_coefficients", &ContactSequence::setFrictionCoefficients,
          py::arg("contact_phase"), py::arg("friction_coefficients"))
    .def("num_impact_events", &ContactSequence::numImpactEvents)
    .def("num_lift_events", &ContactSequence::numLiftEvents)
    .def("num_discrete_events", &ContactSequence::numDiscreteEvents)
    .def("num_contact_phases", &ContactSequence::numContactPhases)
    .def("contact_status", &ContactSequence::contactStatus,
          py::arg("contact_phase"))
    .def("impact_status", &ContactSequence::impactStatus,
          py::arg("impact_index"))
    .def("impact_time", &ContactSequence::impactTime,
          py::arg("impact_index"))
    .def("lift_time", &ContactSequence::liftTime,
          py::arg("lift_index"))
    .def("event_type", &ContactSequence::eventType,
          py::arg("event_index"))
    .def("event_times", &ContactSequence::eventTimes)
    .def("reserve", &ContactSequence::reserve,
         py::arg("reserved_num_discrete_events"))
    .def("reserved_num_discrete_events", &ContactSequence::reservedNumDiscreteEvents)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(ContactSequence)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(ContactSequence);
}

} // namespace python
} // namespace robotoc