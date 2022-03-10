#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/walking_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(walking_foot_step_planner, m) {
  py::class_<WalkingFootStepPlanner, FootStepPlannerBase, 
             std::shared_ptr<WalkingFootStepPlanner>>(m, "WalkingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("biped_robot"))
    .def("set_gait_pattern", &WalkingFootStepPlanner::setGaitPattern,
         py::arg("step_length"), py::arg("yaw_rate"), 
         py::arg("enable_double_support_phase"))
    .def("init", &WalkingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &WalkingFootStepPlanner::plan,
          py::arg("q"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placement", &WalkingFootStepPlanner::contactPlacement,
          py::arg("step"))
    .def("com", &WalkingFootStepPlanner::com,
          py::arg("step"))
    .def("__str__", [](const WalkingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc