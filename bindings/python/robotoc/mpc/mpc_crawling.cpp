#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_crawling.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_crawling, m) {
  py::class_<MPCCrawling>(m, "MPCCrawling")
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("robot"), py::arg("T"), py::arg("N"), py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCCrawling::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCCrawling::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCCrawling::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCCrawling::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCCrawling::getInitialControlInput)
    .def("get_solution", &MPCCrawling::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCCrawling::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCCrawling::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCCrawling::*)() const>(&MPCCrawling::KKTError))
    .def("get_cost_handle", &MPCCrawling::getCostHandle)
    .def("get_config_cost_handle", &MPCCrawling::getConfigCostHandle)
    .def("get_swing_foot_cost_handle", &MPCCrawling::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCCrawling::getCoMCostHandle)
    .def("get_constraints_handle", &MPCCrawling::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCCrawling::getFrictionConeHandle);
}

} // namespace python
} // namespace robotoc