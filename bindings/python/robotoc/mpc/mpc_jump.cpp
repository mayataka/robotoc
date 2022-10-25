#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_jump.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_jump, m) {
  py::class_<MPCJump>(m, "MPCJump")
    .def(py::init<const Robot&, const double, const int>(),
         py::arg("robot"), py::arg("T"), py::arg("N"))
    .def("set_jump_pattern", &MPCJump::setJumpPattern,
         py::arg("foot_step_planner"), py::arg("flying_time"), py::arg("min_flying_time"), 
         py::arg("ground_time"), py::arg("min_ground_time"))
    .def("init", &MPCJump::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("reset", 
          static_cast<void (MPCJump::*)()>(&MPCJump::reset))
    .def("reset", 
          static_cast<void (MPCJump::*)(const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCJump::reset),
          py::arg("q"), py::arg("v"))
    .def("reset", 
          static_cast<void (MPCJump::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&, const SolverOptions&, const bool)>(&MPCJump::reset),
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"), 
          py::arg("sto")=false)
    .def("set_solver_options", &MPCJump::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCJump::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCJump::getInitialControlInput)
    .def("get_solution", &MPCJump::getSolution)
    .def("get_control_policy", &MPCJump::getControlPolicy,
          py::arg("t"))
    .def("KKT_error", 
          static_cast<double (MPCJump::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCJump::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCJump::*)() const>(&MPCJump::KKTError))
    .def("get_cost_handle", &MPCJump::getCostHandle)
    .def("get_config_cost_handle", &MPCJump::getConfigCostHandle)
    .def("get_constraints_handle", &MPCJump::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCJump::getFrictionConeHandle)
    .def("get_solver", &MPCJump::getSolver)
    .def("get_contact_sequence", &MPCJump::getContactSequence)
    .def("set_robot_properties", &MPCJump::setRobotProperties)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(MPCJump);
}

} // namespace python
} // namespace robotoc