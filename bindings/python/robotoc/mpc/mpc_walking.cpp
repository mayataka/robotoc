#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_walking.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_walking, m) {
  py::class_<MPCWalking>(m, "MPCWalking")
    .def(py::init<const OCP&, const int>(),
         py::arg("ocp"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCWalking::setGaitPattern,
         py::arg("vcom"), py::arg("yaw_rate"), py::arg("swing_time"), 
         py::arg("initial_lift_time"))
    .def("init", &MPCWalking::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCWalking::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCWalking::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCWalking::getInitialControlInput)
    .def("KKT_error", 
          static_cast<double (MPCWalking::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCWalking::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCWalking::*)() const>(&MPCWalking::KKTError))
    .def("get_planner", &MPCWalking::getPlanner);
}

} // namespace python
} // namespace robotoc