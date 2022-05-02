#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_flying_trotting.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_flying_trotting, m) {
  py::class_<MPCFlyingTrotting>(m, "MPCFlyingTrotting")
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), 
         py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCFlyingTrotting::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("flying_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCFlyingTrotting::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCFlyingTrotting::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCFlyingTrotting::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCFlyingTrotting::getInitialControlInput)
    .def("get_solution", &MPCFlyingTrotting::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCFlyingTrotting::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCFlyingTrotting::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCFlyingTrotting::*)() const>(&MPCFlyingTrotting::KKTError))
    .def("get_cost_handle", &MPCFlyingTrotting::getCostHandle)
    .def("get_config_cost_handle", &MPCFlyingTrotting::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCFlyingTrotting::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCFlyingTrotting::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCFlyingTrotting::getCoMCostHandle)
    .def("get_constraints_handle", &MPCFlyingTrotting::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCFlyingTrotting::getFrictionConeHandle)
    .def("get_contact_sequence_handle", &MPCFlyingTrotting::getContactSequenceHandle);
}

} // namespace python
} // namespace robotoc