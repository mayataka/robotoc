#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/flying_trot_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(flying_trot_foot_step_planner, m) {
  py::class_<FlyingTrotFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<FlyingTrotFootStepPlanner>>(m, "FlyingTrotFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", 
          static_cast<void (FlyingTrotFootStepPlanner::*)(const Eigen::Vector3d&, const double)>(&FlyingTrotFootStepPlanner::setGaitPattern),
          py::arg("step_length"), py::arg("yaw_step")) 
    .def("set_gait_pattern", 
          static_cast<void (FlyingTrotFootStepPlanner::*)(const Eigen::Vector3d&, const double, const double, const double, const double)>(&FlyingTrotFootStepPlanner::setGaitPattern),
          py::arg("v_com_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("t_swing"), py::arg("t_stance"), py::arg("gain")) 
    .def("init", &FlyingTrotFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &FlyingTrotFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (FlyingTrotFootStepPlanner::*)(const int) const>(&FlyingTrotFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (FlyingTrotFootStepPlanner::*)() const>(&FlyingTrotFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (FlyingTrotFootStepPlanner::*)(const int) const>(&FlyingTrotFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (FlyingTrotFootStepPlanner::*)() const>(&FlyingTrotFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (FlyingTrotFootStepPlanner::*)(const int) const>(&FlyingTrotFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (FlyingTrotFootStepPlanner::*)() const>(&FlyingTrotFootStepPlanner::R))
    .def("__str__", [](const FlyingTrotFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc