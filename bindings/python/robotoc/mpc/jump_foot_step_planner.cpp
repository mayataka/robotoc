#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/jump_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(jump_foot_step_planner, m) {
  py::class_<JumpFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<JumpFootStepPlanner>>(m, "JumpFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    .def("set_jump_pattern", &JumpFootStepPlanner::setJumpPattern,
         py::arg("jump_length"), py::arg("jump_yaw"))
    .def("init", &JumpFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &JumpFootStepPlanner::plan,
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (JumpFootStepPlanner::*)(const int) const>(&JumpFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (JumpFootStepPlanner::*)() const>(&JumpFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (JumpFootStepPlanner::*)(const int) const>(&JumpFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (JumpFootStepPlanner::*)() const>(&JumpFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (JumpFootStepPlanner::*)(const int) const>(&JumpFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (JumpFootStepPlanner::*)() const>(&JumpFootStepPlanner::R))
    .def("__str__", [](const JumpFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc