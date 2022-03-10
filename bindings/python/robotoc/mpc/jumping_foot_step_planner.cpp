#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/jumping_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(jumping_foot_step_planner, m) {
  py::class_<JumpingFootStepPlanner, FootStepPlannerBase, 
             std::shared_ptr<JumpingFootStepPlanner>>(m, "JumpingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("robot"))
    .def("set_jump_pattern", &JumpingFootStepPlanner::setJumpPattern,
         py::arg("jump_length"), py::arg("jump_yaw"))
    .def("init", &JumpingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &JumpingFootStepPlanner::plan,
          py::arg("q"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placement", &JumpingFootStepPlanner::contactPlacement,
          py::arg("step"))
    .def("com", &JumpingFootStepPlanner::com,
          py::arg("step"))
    .def("__str__", [](const JumpingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc