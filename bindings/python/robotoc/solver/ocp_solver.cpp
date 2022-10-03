#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/utils/pybind11_macros.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(ocp_solver, m) {
  py::class_<OCPSolver>(m, "OCPSolver")
    .def(py::init<const OCP&, const SolverOptions&, const int>(),
          py::arg("ocp"), py::arg("solver_options")=SolverOptions(), 
          py::arg("nthreads")=1)
    .def("set_solver_options", &OCPSolver::setSolverOptions,
          py::arg("solver_options"))
    .def("discretize", &OCPSolver::discretize,
          py::arg("t"))
    .def("init_constraints", &OCPSolver::initConstraints)
    .def("update_solution", &OCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("solve", &OCPSolver::solve,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("init_solver")=true)
    .def("get_solver_statistics", &OCPSolver::getSolverStatistics)
    .def("get_solution", 
          static_cast<const Solution& (OCPSolver::*)() const>(&OCPSolver::getSolution))
    .def("get_solution", 
          static_cast<const SplitSolution& (OCPSolver::*)(const int) const>(&OCPSolver::getSolution),
          py::arg("stage"))
    .def("get_solution", 
          static_cast<std::vector<Eigen::VectorXd> (OCPSolver::*)(const std::string&, const std::string&) const>(&OCPSolver::getSolution),
          py::arg("name"), py::arg("option")="")
    .def("get_LQR_policy", &OCPSolver::getLQRPolicy)
    .def("get_riccati_factorization", &OCPSolver::getRiccatiFactorization)
    .def("set_solution", 
          static_cast<void (OCPSolver::*)(const Solution&)>(&OCPSolver::setSolution),
          py::arg("s"))
    .def("set_solution", 
          static_cast<void (OCPSolver::*)(const std::string&, const Eigen::VectorXd&)>(&OCPSolver::setSolution),
          py::arg("name"), py::arg("value"))
    .def("KKT_error", 
          static_cast<double (OCPSolver::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&OCPSolver::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (OCPSolver::*)() const>(&OCPSolver::KKTError))
    .def("get_time_discretization", &OCPSolver::getTimeDiscretization)
    .def("set_robot_properties", &OCPSolver::setRobotProperties,
          py::arg("properties"))
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(OCPSolver)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(OCPSolver);
}

} // namespace python
} // namespace robotoc