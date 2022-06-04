#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/mpc/mpc_pace.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(mpc_pace, m) {
  py::class_<MPCPace>(m, "MPCPace")
    .def(py::init<const Robot&, const double, const int, const int,  const int>(),
         py::arg("quadruped_robot"), py::arg("T"), py::arg("N"), 
         py::arg("max_steps"), py::arg("nthreads"))
    .def("set_gait_pattern", &MPCPace::setGaitPattern,
         py::arg("planner"), py::arg("swing_height"), py::arg("swing_time"), 
         py::arg("stance_time"), py::arg("swing_start_time"))
    .def("init", &MPCPace::init,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("solver_options"))
    .def("set_solver_options", &MPCPace::setSolverOptions,
          py::arg("solver_options"))
    .def("update_solution", &MPCPace::updateSolution,
          py::arg("t"), py::arg("dt"), py::arg("q"), py::arg("v"))
    .def("get_initial_control_input", &MPCPace::getInitialControlInput)
    .def("get_solution", &MPCPace::getSolution)
    .def("KKT_error", 
          static_cast<double (MPCPace::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&MPCPace::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (MPCPace::*)() const>(&MPCPace::KKTError))
    .def("get_cost_handle", &MPCPace::getCostHandle)
    .def("get_config_cost_handle", &MPCPace::getConfigCostHandle)
    .def("get_base_rotation_cost_handle", &MPCPace::getBaseRotationCostHandle)
    .def("get_swing_foot_cost_handle", &MPCPace::getSwingFootCostHandle)
    .def("get_com_cost_handle", &MPCPace::getCoMCostHandle)
    .def("get_constraints_handle", &MPCPace::getConstraintsHandle)
    .def("get_friction_cone_handle", &MPCPace::getFrictionConeHandle)
    .def("get_solver", &MPCPace::getSolver)
    .def("get_contact_sequence", &MPCPace::getContactSequence);
}

} // namespace python
} // namespace robotoc