#include <pybind11/pybind11.h>

#include "robotoc/solver/solver_options.hpp"


namespace robotoc {
namespace python {

namespace py = pybind11;

PYBIND11_MODULE(solver_options, m) {
  py::class_<SolverOptions>(m, "SolverOptions")
    .def(py::init(&SolverOptions::defaultOptions))
    .def_readwrite("max_iter", &SolverOptions::max_iter)
    .def_readwrite("kkt_tol", &SolverOptions::kkt_tol)
    .def_readwrite("mu_init", &SolverOptions::mu_init)
    .def_readwrite("mu_min", &SolverOptions::mu_min)
    .def_readwrite("kkt_tol_mu", &SolverOptions::kkt_tol_mu)
    .def_readwrite("mu_linear_decrease_factor", &SolverOptions::mu_linear_decrease_factor)
    .def_readwrite("mu_superlinear_decrease_power", &SolverOptions::mu_superlinear_decrease_power)
    .def_readwrite("enable_line_search", &SolverOptions::enable_line_search)
    .def_readwrite("line_search_settings", &SolverOptions::line_search_settings)
    .def_readwrite("discretization_method", &SolverOptions::discretization_method)
    .def_readwrite("initial_sto_reg_iter", &SolverOptions::initial_sto_reg_iter)
    .def_readwrite("initial_sto_reg", &SolverOptions::initial_sto_reg)
    .def_readwrite("kkt_tol_mesh", &SolverOptions::kkt_tol_mesh)
    .def_readwrite("max_dt_mesh", &SolverOptions::max_dt_mesh)
    .def_readwrite("max_dts_riccati", &SolverOptions::max_dts_riccati)
    .def_readwrite("print_level", &SolverOptions::print_level)
    .def("__str__", [](const SolverOptions& self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
      });
}

} // namespace python
} // namespace robotoc