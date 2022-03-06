#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_quadrupedal_jumping.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_quadrupedal_jumping, m) {
  py::class_<MPCQuadrupedalJumping>(m, "MPCQuadrupedalJumping")
    .def(py::init<const OCP&, const int>(),
         py::arg("ocp"), py::arg("nthreads"))
    .def("set_jump_pattern", &MPCQuadrupedalJumping::setJumpPattern,
         py::arg("jump_length"), py::arg("jump_yaw"), 
         py::arg("flying_time"), py::arg("min_flying_time"), 
         py::arg("ground_time"), py::arg("min_ground_time"))
    .def("init", &MPCQuadrupedalJumping::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("reset", &MPCQuadrupedalJumping::reset,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("set_solver_options", &MPCQuadrupedalJumping::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCQuadrupedalJumping::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCQuadrupedalJumping::getInitialControlInput)
    .def("KKT_error", 
          static_cast<double (MPCQuadrupedalJumping::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCQuadrupedalJumping::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCQuadrupedalJumping::*)() const>(&MPCQuadrupedalJumping::KKTError));
}

} // namespace python
} // namespace robotoc