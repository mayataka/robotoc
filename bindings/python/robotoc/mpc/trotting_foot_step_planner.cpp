#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/trotting_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trotting_foot_step_planner, m) {
  py::class_<TrottingFootStepPlanner, FootStepPlannerBase, 
             std::shared_ptr<TrottingFootStepPlanner>>(m, "TrottingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", &TrottingFootStepPlanner::setGaitPattern,
         py::arg("step_length"), py::arg("yaw_rate"), 
         py::arg("enable_stance_phase"))
    .def("init", &TrottingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &TrottingFootStepPlanner::plan,
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (TrottingFootStepPlanner::*)(const int) const>(&TrottingFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (TrottingFootStepPlanner::*)() const>(&TrottingFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (TrottingFootStepPlanner::*)(const int) const>(&TrottingFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (TrottingFootStepPlanner::*)() const>(&TrottingFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (TrottingFootStepPlanner::*)(const int) const>(&TrottingFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (TrottingFootStepPlanner::*)() const>(&TrottingFootStepPlanner::R))
    .def("__str__", [](const TrottingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc