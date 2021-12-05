#include "robotoc/solver/solver_options.hpp"


namespace robotoc {

SolverOptions::SolverOptions(const int _max_iter, const double _kkt_tol, 
                             const double _mu_init, const double _mu_min, 
                             const double _kkt_tol_mu, 
                             const double _mu_linear_decrease_factor, 
                             const double _mu_superlinear_decrease_power, 
                             const bool _enable_line_search, 
                             const LineSearchSettings& _line_search_settings, 
                             const DiscretizationMethod& _discretization_method, 
                             const int _initial_sto_reg_iter, 
                             const double _initial_sto_reg, 
                             const double _kkt_tol_mesh, 
                             const double _max_dt_mesh, 
                             const double _max_dts_riccati, 
                             const int _print_level)
  : max_iter(_max_iter),
    kkt_tol(_kkt_tol),
    mu_init(_mu_init),
    mu_min(_mu_min),
    kkt_tol_mu(_kkt_tol_mu),
    mu_linear_decrease_factor(_mu_linear_decrease_factor),
    mu_superlinear_decrease_power(_mu_superlinear_decrease_power),
    enable_line_search(_enable_line_search),
    line_search_settings(_line_search_settings),
    discretization_method(_discretization_method),
    initial_sto_reg_iter(_initial_sto_reg_iter),
    initial_sto_reg(_initial_sto_reg),
    kkt_tol_mesh(_kkt_tol_mesh),
    max_dt_mesh(_max_dt_mesh),
    max_dts_riccati(_max_dts_riccati),
    print_level(_print_level) {
}


SolverOptions::SolverOptions()
  : max_iter(0),
    kkt_tol(0),
    mu_init(0),
    mu_min(0),
    kkt_tol_mu(0),
    mu_linear_decrease_factor(0),
    mu_superlinear_decrease_power(0),
    enable_line_search(false),
    line_search_settings(),
    discretization_method(),
    initial_sto_reg_iter(),
    initial_sto_reg(),
    kkt_tol_mesh(0),
    max_dt_mesh(0),
    max_dts_riccati(0),
    print_level(0) {
}


SolverOptions::~SolverOptions() {
}


SolverOptions SolverOptions::defaultOptions() {
  SolverOptions options;
  options.max_iter = 100;
  options.kkt_tol = 1.0e-07;
  options.mu_init = 1.0e-03;
  options.mu_min = 1.0e-03;
  options.kkt_tol_mu = 1.0e-07;
  options.mu_linear_decrease_factor = 0.2;
  options.mu_superlinear_decrease_power = 1.5;
  options.enable_line_search = false;
  options.line_search_settings = LineSearchSettings::defaultSettings();
  options.discretization_method = DiscretizationMethod::GridBased;
  options.initial_sto_reg_iter = 0;
  options.initial_sto_reg = 1.0e30;
  options.kkt_tol_mesh = 0.1;
  options.max_dt_mesh = 0;
  options.max_dts_riccati = 0.1;
  options.print_level = 1;
  return options;
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
  os << "  mex_dts_riccati: " << max_dts_riccati << std::endl;
  os << "  print_level: " << print_level << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SolverOptions& solver_options) {
  solver_options.disp(os);
  return os;
}

} // namespace robotoc