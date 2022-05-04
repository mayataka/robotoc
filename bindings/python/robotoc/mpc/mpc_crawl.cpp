#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_crawl.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_crawl, m) {
  py::class_<MPCCrawl>(m, "MPCCrawl")
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("robot"), py::arg("T"), py::arg("N"), py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCCrawl::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCCrawl::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCCrawl::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCCrawl::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCCrawl::getInitialControlInput)
    .def("get_solution", &MPCCrawl::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCCrawl::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCCrawl::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCCrawl::*)() const>(&MPCCrawl::KKTError))
    .def("get_cost_handle", &MPCCrawl::getCostHandle)
    .def("get_config_cost_handle", &MPCCrawl::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCCrawl::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCCrawl::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCCrawl::getCoMCostHandle)
    .def("get_constraints_handle", &MPCCrawl::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCCrawl::getFrictionConeHandle);
}

} // namespace python
} // namespace robotoc