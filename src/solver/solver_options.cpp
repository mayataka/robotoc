#include "robotoc/solver/solver_options.hpp"


namespace robotoc {

void SolverOptions::disp(std::ostream& os) const {
  os << "Solver options:" << "\n";
  os << "  nthreads: " << nthreads << "\n";
  os << "  max_iter: " << max_iter << "\n";
  os << "  kkt_tol: " << kkt_tol << "\n";
  os << "  mu_init: " << mu_init << "\n";
  os << "  mu_min: " << mu_min << "\n";
  os << "  kkt_tol_mu: " << kkt_tol_mu << "\n";
  os << "  mu_linear_decrease_factor: " << mu_linear_decrease_factor << "\n";
  os << "  mu_superlinear_decrease_power: " << mu_superlinear_decrease_power << "\n";
  os << "  enable_line_search: " << std::boolalpha << enable_line_search << "\n";
  os << "  line_search_settings: " << line_search_settings << "\n";
  os << "  discretization_method: ";
  if (discretization_method == DiscretizationMethod::GridBased) os << "GridBased" << "\n";
  else os << "PhaseBased" << "\n";
  os << "  initial_sto_reg_iter: " << initial_sto_reg_iter << "\n";
  os << "  initial_sto_reg: " << initial_sto_reg << "\n";
  os << "  kkt_tol_mesh: " << kkt_tol_mesh << "\n";
  os << "  max_dt_mesh: " << max_dt_mesh << "\n";
  os << "  mex_dts_riccati: " << max_dts_riccati << "\n";
  os << "  enable_solution_interpolation: " << std::boolalpha << enable_solution_interpolation << "\n";
  os << "  interpolation_order: ";
  if (interpolation_order == InterpolationOrder::Linear) os << "Linear" << "\n";
  else os << "Zero" << "\n";
  os << "  enable_benchmark: " << std::boolalpha << enable_benchmark << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverOptions& solver_options) {
  solver_options.disp(os);
  return os;
}

} // namespace robotoc