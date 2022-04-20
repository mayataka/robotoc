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
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), 
         py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCWalking::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("double_support_time"), py::arg("swing_start_time"))
    .def("init", &MPCWalking::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCWalking::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCWalking::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCWalking::getInitialControlInput)
    .def("get_solution", &MPCWalking::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCWalking::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCWalking::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCWalking::*)() const>(&MPCWalking::KKTError))
    .def("get_cost_handle", &MPCWalking::getCostHandle)
    .def("get_config_cost_handle", &MPCWalking::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCWalking::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCWalking::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCWalking::getCoMCostHandle)
    .def("get_constraints_handle", &MPCWalking::getConstraintsHandle)
    .def("get_wrench_cone_handle", &MPCWalking::getWrenchConeHandle)
    .def("get_impulse_wrench_cone_handle", &MPCWalking::getImpulseWrenchConeHandle);
}

} // namespace python
} // namespace robotoc