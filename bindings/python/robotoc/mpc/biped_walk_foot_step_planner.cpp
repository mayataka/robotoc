#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/biped_walk_foot_step_planner.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


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
    .def("contact_placements", &BipedWalkFootStepPlanner::contactPlacements,
          py::arg("step"))
    .def("com", &BipedWalkFootStepPlanner::CoM,
          py::arg("step"))
    .def("R", &BipedWalkFootStepPlanner::R,
          py::arg("step"))
    .def("size", &BipedWalkFootStepPlanner::size)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(BipedWalkFootStepPlanner)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(BipedWalkFootStepPlanner);
}

} // namespace python
} // namespace robotoc