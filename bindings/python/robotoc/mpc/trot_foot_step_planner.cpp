#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/trot_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(trot_foot_step_planner, m) {
  py::class_<TrotFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<TrotFootStepPlanner>>(m, "TrotFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", 
          static_cast<void (TrotFootStepPlanner::*)(const Eigen::Vector3d&, const double, const bool)>(&TrotFootStepPlanner::setGaitPattern),
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_stance_phase")) 
    .def("set_gait_pattern", 
          static_cast<void (TrotFootStepPlanner::*)(const Eigen::Vector3d&, const double, const double, const double, const double)>(&TrotFootStepPlanner::setGaitPattern),
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("swing_time"), py::arg("stance_time"), py::arg("gain")) 
    .def("init", &TrotFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &TrotFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (TrotFootStepPlanner::*)(const int) const>(&TrotFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (TrotFootStepPlanner::*)() const>(&TrotFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (TrotFootStepPlanner::*)(const int) const>(&TrotFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (TrotFootStepPlanner::*)() const>(&TrotFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (TrotFootStepPlanner::*)(const int) const>(&TrotFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (TrotFootStepPlanner::*)() const>(&TrotFootStepPlanner::R))
    .def("__str__", [](const TrotFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc