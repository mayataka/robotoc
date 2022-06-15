#include "robotoc/solver/solver_options.hpp"


namespace robotoc {

SolverOptions::SolverOptions() {
  max_iter = 100;
  kkt_tol = 1.0e-07;
  mu_init = 1.0e-03;
  mu_min = 1.0e-03;
  kkt_tol_mu = 1.0e-07;
  mu_linear_decrease_factor = 0.2;
  mu_superlinear_decrease_power = 1.5;
  enable_line_search = false;
  line_search_settings = LineSearchSettings::defaultSettings();
  discretization_method = DiscretizationMethod::GridBased;
  initial_sto_reg_iter = 0;
  initial_sto_reg = 1.0e30;
  kkt_tol_mesh = 0.1;
  max_dt_mesh = 0;
  max_dts_riccati = 0.1;
  enable_benchmark = false;
}


SolverOptions::~SolverOptions() {
}


SolverOptions SolverOptions::defaultOptions() {
  return SolverOptions();
}


void SolverOptions::disp(std::ostream& os) const {
  os << "Solver options:" << std::endl;
  os << "  max_iter: " << max_iter << std::endl;
  os << "  kkt_tol: " << kkt_tol << std::endl;
  os << "  mu_init: " << mu_init << std::endl;
  os << "  mu_min: " << mu_min << std::endl;
  os << "  kkt_tol_mu: " << kkt_tol_mu << std::endl;
  os << "  mu_linear_decrease_factor: " << mu_linear_decrease_factor << std::endl;
  os << "  mu_superlinear_decrease_power: " << mu_superlinear_decrease_power << std::endl;
  os << "  enable_line_search: " << std::boolalpha << enable_line_search << std::endl;
  os << "  line_search_settings: " << line_search_settings << std::endl;
  os << "  discretization_method: ";
  if (discretization_method == DiscretizationMethod::GridBased) os << "grid-based" << std::endl;
  else os << "phase-based" << std::endl;
  os << "  initial_sto_reg_iter: " << initial_sto_reg_iter << std::endl;
  os << "  initial_sto_reg: " << initial_sto_reg << std::endl;
  os << "  kkt_tol_mesh: " << kkt_tol_mesh << std::endl;
  os << "  max_dt_mesh: " << max_dt_mesh << std::endl;
  os << "  mex_dts_riccati: " << max_dts_riccati << std::flush;
  os << "  enable_benchmark: " << std::boolalpha << enable_benchmark << std::endl;
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverOptions& solver_options) {
  solver_options.disp(os);
  return os;
}

} // namespace robotoc