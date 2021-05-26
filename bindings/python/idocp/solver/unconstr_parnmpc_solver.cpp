#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/solver/unconstr_parnmpc_solver.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(unconstr_parnmpc_solver, m) {
  py::class_<UnconstrParNMPCSolver>(m, "UnconstrParNMPCSolver")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int, 
                  const int>(),
         py::arg("robot"), py::arg("cost"), py::arg("constraints"),
         py::arg("T"), py::arg("N"), py::arg("nthreads")=1)
    .def("init_constraints", &UnconstrParNMPCSolver::initConstraints)
    .def("update_solution", &UnconstrParNMPCSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), 
          py::arg("line_search")=false)
    .def("get_solution", static_cast<const SplitSolution& (UnconstrParNMPCSolver::*)(const int stage) const>(&UnconstrParNMPCSolver::getSolution))
    .def("get_solution", static_cast<std::vector<Eigen::VectorXd> (UnconstrParNMPCSolver::*)(const std::string&) const>(&UnconstrParNMPCSolver::getSolution))
    .def("set_solution", &UnconstrParNMPCSolver::setSolution)
    .def("compute_KKT_residual", &UnconstrParNMPCSolver::computeKKTResidual)
    .def("KKT_error", &UnconstrParNMPCSolver::KKTError);
}

} // namespace python
} // namespace idocp