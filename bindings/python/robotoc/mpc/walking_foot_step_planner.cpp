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
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (WalkingFootStepPlanner::*)(const int) const>(&WalkingFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (WalkingFootStepPlanner::*)() const>(&WalkingFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (WalkingFootStepPlanner::*)(const int) const>(&WalkingFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (WalkingFootStepPlanner::*)() const>(&WalkingFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (WalkingFootStepPlanner::*)(const int) const>(&WalkingFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (WalkingFootStepPlanner::*)() const>(&WalkingFootStepPlanner::R))
    .def("__str__", [](const WalkingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc