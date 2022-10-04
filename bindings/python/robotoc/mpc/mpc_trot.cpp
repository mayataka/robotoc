#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_trot.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_trot, m) {
  py::class_<MPCTrot>(m, "MPCTrot")
    .def(py::init<const Robot&, const double, const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"))
    .def("set_gait_pattern", &MPCTrot::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCTrot::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("reset", 
          static_cast<void (MPCTrot::*)()>(&MPCTrot::reset))
    .def("reset", 
          static_cast<void (MPCTrot::*)(const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCTrot::reset),
          py::arg("q"), py::arg("v"))
    .def("set_solver_options", &MPCTrot::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCTrot::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCTrot::getInitialControlInput)
    .def("get_solution", &MPCTrot::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCTrot::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCTrot::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCTrot::*)() const>(&MPCTrot::KKTError))
    .def("get_cost_handle", &MPCTrot::getCostHandle)
    .def("get_config_cost_handle", &MPCTrot::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCTrot::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCTrot::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCTrot::getCoMCostHandle)
    .def("get_constraints_handle", &MPCTrot::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCTrot::getFrictionConeHandle)
    .def("get_solver", &MPCTrot::getSolver)
    .def("get_contact_sequence", &MPCTrot::getContactSequence)
    .def("set_robot_properties", &MPCTrot::setRobotProperties)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(MPCTrot);
}

} // namespace python
} // namespace robotoc