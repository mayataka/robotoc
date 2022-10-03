#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/solver/unconstr_ocp_solver.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_ocp_solver, m) {
  py::class_<UnconstrOCPSolver>(m, "UnconstrOCPSolver")
    .def(py::init<const OCP&, const SolverOptions&, const int>(),
          py::arg("ocp"), py::arg("solver_options")=SolverOptions(), 
          py::arg("nthreads")=1)
    .def("set_solver_options", &UnconstrOCPSolver::setSolverOptions,
          py::arg("solver_options"))
    .def("init_constraints", &UnconstrOCPSolver::initConstraints)
    .def("update_solution", &UnconstrOCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("solve", &UnconstrOCPSolver::solve,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("init_solver")=true)
    .def("get_solver_statistics", &UnconstrOCPSolver::getSolverStatistics)
    .def("get_solution", 
          static_cast<const SplitSolution& (UnconstrOCPSolver::*)(const int) const>(&UnconstrOCPSolver::getSolution))
    .def("get_solution", 
          static_cast<std::vector<Eigen::VectorXd> (UnconstrOCPSolver::*)(const std::string&) const>(&UnconstrOCPSolver::getSolution))
    .def("set_solution", &UnconstrOCPSolver::setSolution,
          py::arg("name"), py::arg("value"))
    .def("get_LQR_policy", &UnconstrOCPSolver::getLQRPolicy)
    .def("KKT_error", 
          static_cast<double (UnconstrOCPSolver::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&UnconstrOCPSolver::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (UnconstrOCPSolver::*)() const>(&UnconstrOCPSolver::KKTError))
    .def("get_time_discretization", &UnconstrOCPSolver::getTimeDiscretization)
    .def("set_robot_properties", &UnconstrOCPSolver::setRobotProperties,
          py::arg("properties"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(UnconstrOCPSolver);
}

} // namespace python
} // namespace robotoc