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
    .def("set_gait_pattern", &CrawlFootStepPlanner::setGaitPattern,
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_stance_phase")) 
    .def("set_raibert_gait_pattern", &CrawlFootStepPlanner::setRaibertGaitPattern,
          py::arg("vcom_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("swing_time"), py::arg("stance_time"), py::arg("gain")) 
    .def("init", &CrawlFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &CrawlFootStepPlanner::plan,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (CrawlFootStepPlanner::*)(const int) const>(&CrawlFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (CrawlFootStepPlanner::*)() const>(&CrawlFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (CrawlFootStepPlanner::*)(const int) const>(&CrawlFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (CrawlFootStepPlanner::*)() const>(&CrawlFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (CrawlFootStepPlanner::*)(const int) const>(&CrawlFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (CrawlFootStepPlanner::*)() const>(&CrawlFootStepPlanner::R))
    .def("__str__", [](const CrawlFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc