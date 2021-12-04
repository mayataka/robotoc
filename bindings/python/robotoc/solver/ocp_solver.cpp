#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/solver/ocp_solver.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(ocp_solver, m) {
  py::class_<OCPSolver>(m, "OCPSolver")
    .def(py::init<const OCP&, const std::shared_ptr<ContactSequence>&, 
                  const SolverOptions&, const int>(),
          py::arg("ocp"), py::arg("contact_sequence"), 
          py::arg("solver_options")=SolverOptions::defaultOptions(), 
          py::arg("nthreads")=1)
    .def("set_solver_options", &OCPSolver::setSolverOptions,
          py::arg("solver_options"))
    .def("mesh_refinement", &OCPSolver::meshRefinement,
          py::arg("t"))
    .def("init_constraints", &OCPSolver::initConstraints,
          py::arg("t"))
    .def("update_solution", &OCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("solve", &OCPSolver::solve,
          py::arg("t"), py::arg("q"), py::arg("v"), py::arg("init_solver")=false)
    .def("get_solution", 
          static_cast<const SplitSolution& (OCPSolver::*)(const int) const>(&OCPSolver::getSolution))
    .def("get_solution", 
          static_cast<std::vector<Eigen::VectorXd> (OCPSolver::*)(const std::string&, const std::string&) const>(&OCPSolver::getSolution),
          py::arg("name"), py::arg("option")="")
    .def("set_solution", &OCPSolver::setSolution,
          py::arg("name"), py::arg("value"))
    .def("KKT_error", 
          static_cast<double (OCPSolver::*)(const double, const Eigen::VectorXd&, const Eigen::VectorXd&)>(&OCPSolver::KKTError),
          py::arg("t"), py::arg("q"), py::arg("v"))
    .def("KKT_error", 
          static_cast<double (OCPSolver::*)() const>(&OCPSolver::KKTError))
    .def("cost", &OCPSolver::cost)
    .def("is_current_solution_feasible", &OCPSolver::isCurrentSolutionFeasible,
          py::arg("verbose")=false)
    .def("get_OCP_discretization", &OCPSolver::getOCPDiscretization)
    .def("set_line_search_settings", &OCPSolver::setLineSearchSettings)
    .def("__str__", [](const OCPSolver& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc