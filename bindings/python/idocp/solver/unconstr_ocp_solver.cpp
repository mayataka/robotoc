#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/solver/unconstr_ocp_solver.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_ocp_solver, m) {
  py::class_<UnconstrOCPSolver>(m, "UnconstrOCPSolver")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int, 
                  const int>(),
         py::arg("robot"), py::arg("cost"), py::arg("constraints"),
         py::arg("T"), py::arg("N"), py::arg("nthreads")=1)
    .def("init_constraints", &UnconstrOCPSolver::initConstraints)
    .def("update_solution", &UnconstrOCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), 
          py::arg("line_search")=false)
    .def("get_solution", static_cast<const SplitSolution& (UnconstrOCPSolver::*)(const int stage) const>(&UnconstrOCPSolver::getSolution))
    .def("get_solution", static_cast<std::vector<Eigen::VectorXd> (UnconstrOCPSolver::*)(const std::string&) const>(&UnconstrOCPSolver::getSolution))
    .def("set_solution", &UnconstrOCPSolver::setSolution)
    .def("compute_KKT_residual", &UnconstrOCPSolver::computeKKTResidual)
    .def("KKT_error", &UnconstrOCPSolver::KKTError)
    .def("cost", &UnconstrOCPSolver::cost);
}

} // namespace python
} // namespace idocp