#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/raibert_flying_trotting_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(raibert_flying_trotting_foot_step_planner, m) {
  py::class_<RaibertFlyingTrottingFootStepPlanner, FootStepPlannerBase, 
             std::shared_ptr<RaibertFlyingTrottingFootStepPlanner>>(m, "RaibertFlyingTrottingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", &RaibertFlyingTrottingFootStepPlanner::setGaitPattern,
         py::arg("v_com_cmd"), py::arg("yaw_rate_cmd"), 
         py::arg("t_swing"), py::arg("t_stance"), py::arg("gain"))
    .def("init", &RaibertFlyingTrottingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &RaibertFlyingTrottingFootStepPlanner::plan,
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (RaibertFlyingTrottingFootStepPlanner::*)(const int) const>(&RaibertFlyingTrottingFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (RaibertFlyingTrottingFootStepPlanner::*)() const>(&RaibertFlyingTrottingFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (RaibertFlyingTrottingFootStepPlanner::*)(const int) const>(&RaibertFlyingTrottingFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (RaibertFlyingTrottingFootStepPlanner::*)() const>(&RaibertFlyingTrottingFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (RaibertFlyingTrottingFootStepPlanner::*)(const int) const>(&RaibertFlyingTrottingFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (RaibertFlyingTrottingFootStepPlanner::*)() const>(&RaibertFlyingTrottingFootStepPlanner::R))
    .def("__str__", [](const RaibertFlyingTrottingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc