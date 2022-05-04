#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/crawling_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(crawling_foot_step_planner, m) {
  py::class_<CrawlingFootStepPlanner, ContactPlannerBase, 
             std::shared_ptr<CrawlingFootStepPlanner>>(m, "CrawlingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", 
          static_cast<void (CrawlingFootStepPlanner::*)(const Eigen::Vector3d&, const double, const bool)>(&CrawlingFootStepPlanner::setGaitPattern),
          py::arg("step_length"), py::arg("step_yaw"), py::arg("enable_stance_phase")) 
    .def("set_gait_pattern", 
          static_cast<void (CrawlingFootStepPlanner::*)(const Eigen::Vector3d&, const double, const double, const double, const double)>(&CrawlingFootStepPlanner::setGaitPattern),
          py::arg("v_com_cmd"), py::arg("yaw_rate_cmd"), 
          py::arg("t_swing"), py::arg("t_stance"), py::arg("gain")) 
    .def("init", &CrawlingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &CrawlingFootStepPlanner::plan,
          py::arg("q"), py::arg("v"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", 
          static_cast<const std::vector<Eigen::Vector3d>& (CrawlingFootStepPlanner::*)(const int) const>(&CrawlingFootStepPlanner::contactPosition),
          py::arg("step"))
    .def("contact_position", 
          static_cast<const std::vector<std::vector<Eigen::Vector3d>>& (CrawlingFootStepPlanner::*)() const>(&CrawlingFootStepPlanner::contactPosition))
    .def("com", 
          static_cast<const Eigen::Vector3d& (CrawlingFootStepPlanner::*)(const int) const>(&CrawlingFootStepPlanner::com),
          py::arg("step"))
    .def("com", 
          static_cast<const std::vector<Eigen::Vector3d>& (CrawlingFootStepPlanner::*)() const>(&CrawlingFootStepPlanner::com))
    .def("R", 
          static_cast<const Eigen::Matrix3d& (CrawlingFootStepPlanner::*)(const int) const>(&CrawlingFootStepPlanner::R),
          py::arg("step"))
    .def("R", 
          static_cast<const std::vector<Eigen::Matrix3d>& (CrawlingFootStepPlanner::*)() const>(&CrawlingFootStepPlanner::R))
    .def("__str__", [](const CrawlingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc