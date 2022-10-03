#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_biped_walk.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_biped_walk, m) {
  py::class_<MPCBipedWalk>(m, "MPCBipedWalk")
    .def(py::init<const Robot&, const double, const int,  const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCBipedWalk::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("double_support_time"), py::arg("swing_start_time"))
    .def("init", &MPCBipedWalk::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("reset", 
          static_cast<void (MPCBipedWalk::*)()>(&MPCBipedWalk::reset))
    .def("reset", 
          static_cast<void (MPCBipedWalk::*)(const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCBipedWalk::reset),
          py::arg("q"), py::arg("v"))
    .def("set_solver_options", &MPCBipedWalk::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCBipedWalk::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCBipedWalk::getInitialControlInput)
    .def("get_solution", &MPCBipedWalk::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCBipedWalk::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCBipedWalk::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCBipedWalk::*)() const>(&MPCBipedWalk::KKTError))
    .def("get_cost_handle", &MPCBipedWalk::getCostHandle)
    .def("get_config_cost_handle", &MPCBipedWalk::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCBipedWalk::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCBipedWalk::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCBipedWalk::getCoMCostHandle)
    .def("get_constraints_handle", &MPCBipedWalk::getConstraintsHandle)
    .def("get_contact_wrench_cone_handle", &MPCBipedWalk::getContactWrenchConeHandle)
    .def("get_impulse_wrench_cone_handle", &MPCBipedWalk::getImpulseWrenchConeHandle)
    .def("get_solver", &MPCBipedWalk::getSolver)
    .def("get_contact_sequence", &MPCBipedWalk::getContactSequence)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(MPCBipedWalk);
}

} // namespace python
} // namespace robotoc