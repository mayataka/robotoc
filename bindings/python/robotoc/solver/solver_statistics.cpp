#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "robotoc/solver/solver_statistics.hpp"
#include "robotoc/utils/pybind11_macros.hpp"



namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(solver_statistics, m) {
  py::class_<SolverStatistics>(m, "SolverStatistics")
    .def(py::init<>())
    .def_readonly("convergence", &SolverStatistics::convergence)
    .def_readonly("iter", &SolverStatistics::iter)
    .def_readonly("kkt_error", &SolverStatistics::kkt_error)
    .def_readonly("primal_step_size", &SolverStatistics::primal_step_size)
    .def_readonly("dual_step_size", &SolverStatistics::dual_step_size)
    .def_readonly("ts", &SolverStatistics::ts)
    .def_readonly("mesh_refinement_iter", &SolverStatistics::mesh_refinement_iter)
    .def_readonly("cpu_time", &SolverStatistics::cpu_time)
    DEFINE_ROBOTOC_PYBIND11_CLASS_CLONE(SolverStatistics)
    DEFINE_ROBOTOC_PYBIND11_CLASS_PRINT(SolverStatistics);
}

} // namespace python
} // namespace robotoc