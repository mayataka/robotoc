#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/biped_walk_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(biped_walk_foot_step_planner, m) {
  py::class_<BipedWalkFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<BipedWalkFootStepPlanner>>(m, "BipedWalkFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("biped_robot"))
    .def("set_gait_pattern", &BipedWalkFootStepPlanner::setGaitPattern,
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_double_support_phase")) 
    .def("set_raibert_gait_pattern", &BipedWalkFootStepPlanner::setRaibertGaitPattern,
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("swing_time"), py::arg("double_support_time"), py::arg("gain")) 
    .def("init", &BipedWalkFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &BipedWalkFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placements", 
          static_cast<const aligned_vector<SE3>& (BipedWalkFootStepPlanner::*)(const int) const>(&BipedWalkFootStepPlanner::contactPlacements),
          py::arg("step"))
    .def("contact_placements", 
          static_cast<const aligned_vector<aligned_vector<SE3>>& (BipedWalkFootStepPlanner::*)() const>(&BipedWalkFootStepPlanner::contactPlacements))
    .def("com", 
          static_cast<const Eigen::Vector3d& (BipedWalkFootStepPlanner::*)(const int) const>(&BipedWalkFootStepPlanner::CoM),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (BipedWalkFootStepPlanner::*)() const>(&BipedWalkFootStepPlanner::CoM))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (BipedWalkFootStepPlanner::*)(const int) const>(&BipedWalkFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (BipedWalkFootStepPlanner::*)() const>(&BipedWalkFootStepPlanner::R))
    .def("__str__", [](const BipedWalkFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc