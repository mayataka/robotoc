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
    .def("clone", [](const JumpFootStepPlanner& self) {
       auto other = self;
       return other;
     })
    .def("set_jump_pattern", &JumpFootStepPlanner::setJumpPattern,
         py::arg("jump_length"), py::arg("jump_yaw"))
    .def("set_contact_surfaces", 
          static_cast<void (JumpFootStepPlanner::*)(const std::vector<Eigen::Matrix3d>& contact_surfaces)>(&JumpFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("set_contact_surfaces", 
          static_cast<void (JumpFootStepPlanner::*)(const std::vector<std::vector<Eigen::Matrix3d>>& contact_surfaces)>(&JumpFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("init", &JumpFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &JumpFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placements", &JumpFootStepPlanner::contactPlacements,
          py::arg("step"))
    .def("contact_positions", &JumpFootStepPlanner::contactPositions,
          py::arg("step"))
    .def("contact_surfaces", &JumpFootStepPlanner::contactSurfaces,
          py::arg("step"))
    .def("com", &JumpFootStepPlanner::CoM,
          py::arg("step"))
    .def("R", &JumpFootStepPlanner::R,
          py::arg("step"))
    .def("size", &JumpFootStepPlanner::size)
    .def("__str__", [](const JumpFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc