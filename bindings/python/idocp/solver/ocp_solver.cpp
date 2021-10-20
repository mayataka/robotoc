#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "idocp/solver/ocp_solver.hpp"


namespace idocp {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(ocp_solver, m) {
  py::class_<OCPSolver>(m, "OCPSolver")
    .def(py::init<const Robot&, const std::shared_ptr<CostFunction>&,
                  const std::shared_ptr<Constraints>&, const double, const int, 
                  const int, const int>(),
         py::arg("robot"), py::arg("cost"), py::arg("constraints"),
         py::arg("T"), py::arg("N"), py::arg("max_num_impulse")=0,
         py::arg("nthreads")=1)
    .def("init_constraints", &OCPSolver::initConstraints)
    .def("update_solution", &OCPSolver::updateSolution,
          py::arg("t"), py::arg("q"), py::arg("v"), 
          py::arg("line_search")=false)
    .def("get_solution", static_cast<const SplitSolution& (OCPSolver::*)(const int stage) const>(&OCPSolver::getSolution))
    .def("get_solution", static_cast<std::vector<Eigen::VectorXd> (OCPSolver::*)(const std::string&, const std::string&) const>(&OCPSolver::getSolution),
          py::arg("name"), py::arg("option")="")
    .def("set_solution", &OCPSolver::setSolution)
    .def("set_contact_status_uniformly", &OCPSolver::setContactStatusUniformly)
    .def("set_contact_points", &OCPSolver::setContactPoints)
    .def("push_back_contact_status", &OCPSolver::pushBackContactStatus)
    .def("pop_back_contact_status", &OCPSolver::popBackContactStatus,
          py::arg("t"), py::arg("extrapolate_solution")=false)
    .def("pop_front_contact_status", &OCPSolver::popFrontContactStatus,
          py::arg("t"), py::arg("extrapolate_solution")=false)
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
} // namespace idocp