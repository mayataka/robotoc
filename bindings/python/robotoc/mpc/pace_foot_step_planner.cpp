#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/pace_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(pace_foot_step_planner, m) {
  py::class_<PaceFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<PaceFootStepPlanner>>(m, "PaceFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("clone", &PaceFootStepPlanner::clone)
    .def("set_gait_pattern", &PaceFootStepPlanner::setGaitPattern,
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_stance_phase")) 
    .def("set_raibert_gait_pattern", &PaceFootStepPlanner::setRaibertGaitPattern,
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("swing_time"), py::arg("stance_time"), py::arg("gain")) 
    .def("set_contact_surfaces", 
          static_cast<void (PaceFootStepPlanner::*)(const std::vector<Eigen::Matrix3d>& contact_surfaces)>(&PaceFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("set_contact_surfaces", 
          static_cast<void (PaceFootStepPlanner::*)(const std::vector<std::vector<Eigen::Matrix3d>>& contact_surfaces)>(&PaceFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("init", &PaceFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &PaceFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placements", &PaceFootStepPlanner::contactPlacements,
          py::arg("step"))
    .def("contact_positions", &PaceFootStepPlanner::contactPositions,
          py::arg("step"))
    .def("contact_surfaces", &PaceFootStepPlanner::contactSurfaces,
          py::arg("step"))
    .def("com", &PaceFootStepPlanner::CoM,
          py::arg("step"))
    .def("R", &PaceFootStepPlanner::R,
          py::arg("step"))
    .def("size", &PaceFootStepPlanner::size)
    .def("__str__", [](const PaceFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc