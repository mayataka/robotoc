#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_flying_trot.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_flying_trot, m) {
  py::class_<MPCFlyingTrot>(m, "MPCFlyingTrot")
    .def(py::init<const Robot&, const double, const int, const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), py::arg("nthreads"))
    .def("clone", &MPCFlyingTrot::clone)
    .def("set_gait_pattern", &MPCFlyingTrot::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("flying_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCFlyingTrot::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("reset", 
          static_cast<void (MPCFlyingTrot::*)()>(&MPCFlyingTrot::reset))
    .def("reset", 
          static_cast<void (MPCFlyingTrot::*)(const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCFlyingTrot::reset),
          py::arg("q"), py::arg("v"))
    .def("set_solver_options", &MPCFlyingTrot::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCFlyingTrot::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCFlyingTrot::getInitialControlInput)
    .def("get_solution", &MPCFlyingTrot::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCFlyingTrot::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCFlyingTrot::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCFlyingTrot::*)() const>(&MPCFlyingTrot::KKTError))
    .def("get_cost_handle", &MPCFlyingTrot::getCostHandle)
    .def("get_config_cost_handle", &MPCFlyingTrot::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCFlyingTrot::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCFlyingTrot::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCFlyingTrot::getCoMCostHandle)
    .def("get_constraints_handle", &MPCFlyingTrot::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCFlyingTrot::getFrictionConeHandle)
    .def("get_contact_sequence_handle", &MPCFlyingTrot::getContactSequenceHandle)
    .def("get_solver", &MPCFlyingTrot::getSolver)
    .def("get_contact_sequence", &MPCFlyingTrot::getContactSequence)
    .def("set_robot_properties", &MPCFlyingTrot::setRobotProperties);
}

} // namespace python
} // namespace robotoc