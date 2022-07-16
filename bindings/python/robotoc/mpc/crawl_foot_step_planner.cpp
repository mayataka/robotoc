#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/crawl_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(crawl_foot_step_planner, m) {
  py::class_<CrawlFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<CrawlFootStepPlanner>>(m, "CrawlFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("clone", [](const CrawlFootStepPlanner& self) {
       auto other = self;
       return other;
     })
    .def("set_gait_pattern", &CrawlFootStepPlanner::setGaitPattern,
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_stance_phase")) 
    .def("set_raibert_gait_pattern", &CrawlFootStepPlanner::setRaibertGaitPattern,
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("swing_time"), py::arg("stance_time"), py::arg("gain")) 
    .def("set_contact_surfaces", 
          static_cast<void (CrawlFootStepPlanner::*)(const std::vector<Eigen::Matrix3d>& contact_surfaces)>(&CrawlFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("set_contact_surfaces", 
          static_cast<void (CrawlFootStepPlanner::*)(const std::vector<std::vector<Eigen::Matrix3d>>& contact_surfaces)>(&CrawlFootStepPlanner::setContactSurfaces),
          py::arg("contact_surfaces"))
    .def("init", &CrawlFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &CrawlFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_placements", &CrawlFootStepPlanner::contactPlacements,
          py::arg("step"))
    .def("contact_positions", &CrawlFootStepPlanner::contactPositions,
          py::arg("step"))
    .def("contact_surfaces", &CrawlFootStepPlanner::contactSurfaces,
          py::arg("step"))
    .def("com", &CrawlFootStepPlanner::CoM,
          py::arg("step"))
    .def("R", &CrawlFootStepPlanner::R,
          py::arg("step"))
    .def("size", &CrawlFootStepPlanner::size)
    .def("__str__", [](const CrawlFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc