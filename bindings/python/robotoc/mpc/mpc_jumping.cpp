#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_jumping.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_jumping, m) {
  py::class_<MPCJumping>(m, "MPCJumping")
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("robot"), py::arg("T"), py::arg("N"), py::arg("max_steps"), py::arg("nthreads"))
    .def("set_jump_pattern", &MPCJumping::setJumpPattern,
         py::arg("foot_step_planner"), py::arg("flying_time"), py::arg("min_flying_time"), 
         py::arg("ground_time"), py::arg("min_ground_time"))
    .def("init", &MPCJumping::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("reset", &MPCJumping::reset,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("set_solver_options", &MPCJumping::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCJumping::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCJumping::getInitialControlInput)
    .def("get_solution", &MPCJumping::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCJumping::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCJumping::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCJumping::*)() const>(&MPCJumping::KKTError))
    .def("get_cost_handle", &MPCJumping::getCostHandle)
    .def("get_config_cost_handle", &MPCJumping::getConfigCostHandle)
    .def("get_constraints_handle", &MPCJumping::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCJumping::getFrictionConeHandle);
}

} // namespace python
} // namespace robotoc