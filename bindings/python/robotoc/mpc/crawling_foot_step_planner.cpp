#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/crawling_foot_step_planner.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(crawling_foot_step_planner, m) {
  py::class_<CrawlingFootStepPlanner, std::shared_ptr<CrawlingFootStepPlanner>>(m, "CrawlingFootStepPlanner")
    .def(py::init<const Robot&>(),
         py::arg("quadruped_robot"))
    .def("set_gait_pattern", &CrawlingFootStepPlanner::setGaitPattern,
         py::arg("step_length"), py::arg("yaw_rate"))
    .def("init", &CrawlingFootStepPlanner::init,
          py::arg("q"))
    .def("plan", &CrawlingFootStepPlanner::plan,
          py::arg("q"), py::arg("contact_status"), py::arg("planning_steps"))
    .def("contact_position", &CrawlingFootStepPlanner::contactPosition,
          py::arg("step"))
    .def("com", &CrawlingFootStepPlanner::com,
          py::arg("step"))
    .def("__str__", [](const CrawlingFootStepPlanner& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc