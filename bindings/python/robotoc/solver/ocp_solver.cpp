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
    .def(py::init<const Robot&, const std::shared_ptr<ContactSequence>&, 
                  const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, 
                  const double, const int, const int>(),
         py::arg("robot"), py::arg("contact_sequence"), py::arg("cost"), 
         py::arg("constraints"), py::arg("T"), py::arg("N"), 
         py::arg("nthreads")=1)
    .def("set_discretization_method", &OCPSolver::setDiscretizationMethod)
    .def("mesh_refinement", &OCPSolver::meshRefinement)
    .def("init_constraints", &OCPSolver::initConstraints)
    .def("update_solution", &OCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), 
          py::arg("line_search")=false)
    .def("get_solution", static_cast<const SplitSolution& (OCPSolver::*)(const int stage) const>(&OCPSolver::getSolution))
    .def("get_solution", static_cast<std::vector<Eigen::VectorXd> (OCPSolver::*)(const std::string&, const std::string&) const>(&OCPSolver::getSolution),
          py::arg("name"), py::arg("option")="")
    .def("set_solution", &OCPSolver::setSolution)
    .def("compute_KKT_residual", &OCPSolver::computeKKTResidual)
    .def("KKT_error", &OCPSolver::KKTError)
    .def("cost", &OCPSolver::cost)
    .def("is_formulation_tractable", &OCPSolver::isFormulationTractable)
    .def("__str__", [](const OCPSolver& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc