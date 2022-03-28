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
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), 
         py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCTrotting::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCTrotting::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCTrotting::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCTrotting::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCTrotting::getInitialControlInput)
    .def("KKT_error", 
          static_cast<double (MPCTrotting::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCTrotting::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCTrotting::*)() const>(&MPCTrotting::KKTError))
    .def("get_cost_handle", &MPCTrotting::getCostHandle)
    .def("get_config_cost_handle", &MPCTrotting::getConfigCostHandle)
    .def("get_swing_foot_cost_handle", &MPCTrotting::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCTrotting::getCoMCostHandle)
    .def("get_constraints_handle", &MPCTrotting::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCTrotting::getFrictionConeHandle);
}

} // namespace python
} // namespace robotoc