#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/trotting_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trotting_foot_step_planner, m) {
  py::class_<TrottingFootStepPlanner, std::shared_ptr<TrottingFootStepPlanner>>(m, "TrottingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", &TrottingFootStepPlanner::setGaitPattern,
         py::arg("step_length"), py::arg("yaw_rate"))
    .def("init", &TrottingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &TrottingFootStepPlanner::plan,
          py::arg("q"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", &TrottingFootStepPlanner::contactPosition,
          py::arg("step"))
    .def("com", &TrottingFootStepPlanner::com,
          py::arg("step"))
    .def("__str__", [](const TrottingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc