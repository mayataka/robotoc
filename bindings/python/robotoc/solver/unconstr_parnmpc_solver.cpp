#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/solver/unconstr_parnmpc_solver.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_parnmpc_solver, m) {
  py::class_<UnconstrParNMPCSolver>(m, "UnconstrParNMPCSolver")
    .def(py::init<const UnconstrParNMPC&, const SolverOptions&, const int>(),
          py::arg("parnmpc"), 
          py::arg("solver_options")=SolverOptions::defaultOptions(), 
          py::arg("nthreads")=1)
    .def("set_solver_options", &UnconstrParNMPCSolver::setSolverOptions,
          py::arg("solver_options"))
    .def("init_constraints", &UnconstrParNMPCSolver::initConstraints)
    .def("init_backward_correction", &UnconstrParNMPCSolver::initBackwardCorrection,
          py::arg("t"))
    .def("update_solution", &UnconstrParNMPCSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("solve", &UnconstrParNMPCSolver::solve,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("init_solver")=false)
    .def("get_solver_statistics", &UnconstrParNMPCSolver::getSolverStatistics)
    .def("get_solution", 
          static_cast<const SplitSolution& (UnconstrParNMPCSolver::*)(const int) const>(&UnconstrParNMPCSolver::getSolution))
    .def("get_solution", 
          static_cast<std::vector<Eigen::VectorXd> (UnconstrParNMPCSolver::*)(const std::string&) const>(&UnconstrParNMPCSolver::getSolution))
    .def("set_solution", &UnconstrParNMPCSolver::setSolution,
          py::arg("name"), py::arg("value"))
    .def("KKT_error", 
          static_cast<double (UnconstrParNMPCSolver::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&UnconstrParNMPCSolver::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (UnconstrParNMPCSolver::*)() const>(&UnconstrParNMPCSolver::KKTError))
    .def("cost", &UnconstrParNMPCSolver::cost);
}

} // namespace python
} // namespace robotoc