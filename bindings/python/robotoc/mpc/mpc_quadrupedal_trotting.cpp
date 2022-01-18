#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_quadrupedal_trotting.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_quadrupedal_trotting, m) {
  py::class_<MPCQuadrupedalTrotting>(m, "MPCQuadrupedalTrotting")
    .def(py::init<const OCP&, const int>(),
         py::arg("ocp"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCQuadrupedalTrotting::setGaitPattern,
         py::arg("step_length"), py::arg("step_height"), py::arg("swing_time"), 
         py::arg("initial_lift_time"))
    .def("init", &MPCQuadrupedalTrotting::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCQuadrupedalTrotting::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCQuadrupedalTrotting::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCQuadrupedalTrotting::getInitialControlInput)
    .def("KKT_error", 
          static_cast<double (MPCQuadrupedalTrotting::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCQuadrupedalTrotting::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCQuadrupedalTrotting::*)() const>(&MPCQuadrupedalTrotting::KKTError));
}

} // namespace python
} // namespace robotoc