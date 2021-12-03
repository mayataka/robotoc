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
                  const int>(),
          py::arg("ocp"), py::arg("contact_sequence"), py::arg("nthreads")=1)
    .def("mesh_refinement", &OCPSolver::meshRefinement)
    .def("init_constraints", &OCPSolver::initConstraints)
    .def("update_solution", &OCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), 
          py::arg("line_search")=false)
    .def("get_solution", 
          static_cast<const SplitSolution& (OCPSolver::*)(const int stage) const>(&OCPSolver::getSolution))
    .def("get_solution", 
          static_cast<std::vector<Eigen::VectorXd> (OCPSolver::*)(const std::string&, const std::string&) const>(&OCPSolver::getSolution),
          py::arg("name"), py::arg("option")="")
    .def("set_solution", &OCPSolver::setSolution)
    .def("compute_KKT_residual", &OCPSolver::computeKKTResidual)
    .def("KKT_error", &OCPSolver::KKTError)
    .def("cost", &OCPSolver::cost)
    .def("is_current_solution_feasible", &OCPSolver::isCurrentSolutionFeasible,
          py::arg("verbose")=false)
    .def("get_OCP_discretization", &OCPSolver::getOCPDiscretization)
    .def("set_STO_regularization", &OCPSolver::setSTORegularization)
    .def("set_line_search_settings", &OCPSolver::setLineSearchSettings)
    .def("__str__", [](const OCPSolver& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc