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
    .def("set_gait_pattern", &FlyingTrotFootStepPlanner::setGaitPattern,
          py::arg("step_length"), py::arg("step_yaw")) 
    .def("set_raibert_gait_pattern", &FlyingTrotFootStepPlanner::setRaibertGaitPattern,
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("flying_time"), py::arg("stance_time"), py::arg("gain")) 
    .def("set_contact_surfaces", 
          static_cast<void (FlyingTrotFootStepPlanner::*)(const std::vector<Eigen::Matrix3d>& contact_surfaces)>(&FlyingTrotFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("set_contact_surfaces", 
          static_cast<void (FlyingTrotFootStepPlanner::*)(const std::vector<std::vector<Eigen::Matrix3d>>& contact_surfaces)>(&FlyingTrotFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("init", &FlyingTrotFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &FlyingTrotFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placements", &FlyingTrotFootStepPlanner::contactPlacements,
          py::arg("step"))
    .def("contact_positions", &FlyingTrotFootStepPlanner::contactPositions,
          py::arg("step"))
    .def("contact_surfaces", &FlyingTrotFootStepPlanner::contactSurfaces,
          py::arg("step"))
    .def("com", &FlyingTrotFootStepPlanner::CoM,
          py::arg("step"))
    .def("R", &FlyingTrotFootStepPlanner::R,
          py::arg("step"))
    .def("size", &FlyingTrotFootStepPlanner::size)
    .def("__str__", [](const FlyingTrotFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc