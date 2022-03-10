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
    .def(py::init<const OCP&, const int>(),
         py::arg("ocp"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCCrawling::setGaitPattern,
         py::arg("planner"), py::arg("swing_time"), py::arg("initial_lift_time"))
    .def("init", &MPCCrawling::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCCrawling::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCCrawling::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCCrawling::getInitialControlInput)
    .def("KKT_error", 
          static_cast<double (MPCCrawling::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCCrawling::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCCrawling::*)() const>(&MPCCrawling::KKTError));
}

} // namespace python
} // namespace robotoc