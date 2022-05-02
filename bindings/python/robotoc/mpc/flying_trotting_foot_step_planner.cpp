#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/flying_trotting_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(flying_trotting_foot_step_planner, m) {
  py::class_<FlyingTrottingFootStepPlanner, FootStepPlannerBase, 
             std::shared_ptr<FlyingTrottingFootStepPlanner>>(m, "FlyingTrottingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", &FlyingTrottingFootStepPlanner::setGaitPattern,
         py::arg("step_length"), py::arg("yaw_rate"))
    .def("init", &FlyingTrottingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &FlyingTrottingFootStepPlanner::plan,
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (FlyingTrottingFootStepPlanner::*)(const int) const>(&FlyingTrottingFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (FlyingTrottingFootStepPlanner::*)() const>(&FlyingTrottingFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (FlyingTrottingFootStepPlanner::*)(const int) const>(&FlyingTrottingFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (FlyingTrottingFootStepPlanner::*)() const>(&FlyingTrottingFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (FlyingTrottingFootStepPlanner::*)(const int) const>(&FlyingTrottingFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (FlyingTrottingFootStepPlanner::*)() const>(&FlyingTrottingFootStepPlanner::R))
    .def("__str__", [](const FlyingTrottingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc