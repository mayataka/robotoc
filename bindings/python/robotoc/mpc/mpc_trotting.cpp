#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_trotting.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_trotting, m) {
  py::class_<MPCTrotting>(m, "MPCTrotting")
    .def(py::init<const OCP&, const int>(),
         py::arg("ocp"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCTrotting::setGaitPattern,
         py::arg("vcom_cmd"), py::arg("yaw_cmd"), py::arg("swing_time"), 
         py::arg("initial_lift_time"))
    .def("init", &MPCTrotting::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCTrotting::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCTrotting::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCTrotting::getInitialControlInput)
    .def("setDiscreteTimeCoMRef", &MPCTrotting::setDiscreteTimeCoMRef)
    .def("KKT_error", 
          static_cast<double (MPCTrotting::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCTrotting::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCTrotting::*)() const>(&MPCTrotting::KKTError));
}

} // namespace python
} // namespace robotoc