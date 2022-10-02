#include "robotoc/solver/ocp_solver.hpp"

#include <stdexcept>
#include <cassert>
#include <algorithm>


namespace robotoc {

OCPSolver::OCPSolver(const OCP& ocp, 
                     const SolverOptions& solver_options, const int nthreads)
  : robots_(nthreads, ocp.robot),
    contact_sequence_(ocp.contact_sequence),
    time_discretization_(ocp.T, ocp.N, ocp.reserved_num_discrete_events),
    dms_(ocp, nthreads),
    sto_(ocp),
    riccati_recursion_(ocp, nthreads, solver_options.max_dts_riccati),
    line_search_(ocp),
    ocp_(ocp),
    kkt_matrix_(ocp.N+1+ocp.reserved_num_discrete_events, SplitKKTMatrix(ocp.robot)),
    kkt_residual_(ocp.N+1+ocp.reserved_num_discrete_events, SplitKKTResidual(ocp.robot)),
    s_(ocp.N+1+ocp.reserved_num_discrete_events, SplitSolution(ocp.robot)),
    d_(ocp.N+1+ocp.reserved_num_discrete_events, SplitDirection(ocp.robot)),
    riccati_factorization_(ocp.N+1+ocp.reserved_num_discrete_events+1, SplitRiccatiFactorization(ocp.robot)),
    solution_interpolator_(ocp.robot, ocp.N, ocp.reserved_num_discrete_events),
    solver_options_(solver_options),
    solver_statistics_(),
    timer_() {
  if (!ocp.cost) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.cost should not be nullptr!");
  }
  if (!ocp.constraints) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.constraints should not be nullptr!");
  }
  if (!ocp.contact_sequence) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.contact_sequence should not be nullptr!");
  }
  if (ocp.T <= 0) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.T must be positive!");
  }
  if (ocp.N <= 0) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.N must be positive!");
  }
  if (ocp.reserved_num_discrete_events< 0) {
    throw std::out_of_range("[OCPSolver] invalid argument: ocp.reserved_num_discrete_events must be non-negative!");
  }
  if (nthreads <= 0) {
    throw std::out_of_range("[OCPSolver] invalid argument: nthreads must be positive!");
  }
  for (auto& e : s_)  { ocp.robot.normalizeConfiguration(e.q); }
  if (ocp.sto_cost && ocp.sto_constraints) {
    solver_options_.discretization_method = DiscretizationMethod::PhaseBased;
  }
  time_discretization_.setDiscretizationMethod(solver_options_.discretization_method);
}


OCPSolver::OCPSolver()
  : robots_(),
    contact_sequence_(),
    time_discretization_(),
    dms_(),
    sto_(),
    riccati_recursion_(),
    line_search_(),
    ocp_(),
    kkt_matrix_(),
    kkt_residual_(),
    s_(),
    d_(),
    riccati_factorization_(),
    solution_interpolator_(),
    solver_options_(),
    solver_statistics_(),
    timer_() {
}


void OCPSolver::setSolverOptions(const SolverOptions& solver_options) {
  solver_options_ = solver_options;
  if (ocp_.sto_cost && ocp_.sto_constraints) {
    solver_options_.discretization_method = DiscretizationMethod::PhaseBased;
  }
  time_discretization_.setDiscretizationMethod(solver_options_.discretization_method);
  riccati_recursion_.setRegularization(solver_options_.max_dts_riccati);
}


void OCPSolver::discretize(const double t) {
  time_discretization_.discretizeGrid(contact_sequence_, t);
  if (solver_options_.discretization_method == DiscretizationMethod::PhaseBased) {
    time_discretization_.discretizePhase(contact_sequence_, t);
  }
  resizeData();
}


void OCPSolver::initConstraints(const double t) {
  discretize(t);
  dms_.initConstraints(robots_, time_discretization_, s_);
  sto_.initConstraints(time_discretization_);
}


void OCPSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  if ((solver_options_.discretization_method == DiscretizationMethod::PhaseBased)
       || (ocp_.sto_cost && ocp_.sto_constraints)) {
    time_discretization_.discretizePhase(contact_sequence_, t);
  }
  dms_.evalKKT(robots_, time_discretization_, q, v, s_, kkt_matrix_, kkt_residual_);
  sto_.evalKKT(time_discretization_, kkt_matrix_, kkt_residual_);
  riccati_recursion_.backwardRiccatiRecursion(time_discretization_, 
                                              kkt_matrix_, kkt_residual_, 
                                              riccati_factorization_);
  dms_.computeInitialStateDirection(robots_[0], q, v, s_, d_);
  riccati_recursion_.forwardRiccatiRecursion(time_discretization_, 
                                             kkt_matrix_, kkt_residual_, 
                                             riccati_factorization_, d_);
  dms_.computeStepSizes(time_discretization_, d_);
  sto_.computeStepSizes(time_discretization_, d_);
  double primal_step_size = std::min(dms_.maxPrimalStepSize(), 
                                     sto_.maxPrimalStepSize());
  const double dual_step_size = std::min(dms_.maxDualStepSize(),
                                         sto_.maxDualStepSize());
  if (solver_options_.enable_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(dms_, robots_, 
                                                    time_discretization_, 
                                                    q, v, s_, d_, 
                                                    max_primal_step_size);
  }
  solver_statistics_.primal_step_size.push_back(primal_step_size);
  solver_statistics_.dual_step_size.push_back(dual_step_size);
  dms_.integrateSolution(robots_, time_discretization_, 
                         primal_step_size, dual_step_size, kkt_matrix_, d_, s_);
  sto_.integrateSolution(time_discretization_, primal_step_size, dual_step_size, d_);
} 


void OCPSolver::solve(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, const bool init_solver) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[OCPSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[OCPSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
  }
  if (solver_options_.enable_benchmark) {
    timer_.tick();
  }
  if (init_solver) {
    discretize(t);
    if (solver_options_.enable_solution_interpolation) {
      solution_interpolator_.interpolateTimeBased(robots_[0], time_discretization_, s_);
    }
    dms_.initConstraints(robots_, time_discretization_, s_);
    sto_.initConstraints(time_discretization_);
    line_search_.clearFilter();
  }
  solver_statistics_.clear(); 
  int inner_iter = 0;
  for (int iter=0; iter<solver_options_.max_iter; ++iter, ++inner_iter) {
    if (ocp_.sto_cost && ocp_.sto_constraints) {
      if (inner_iter < solver_options_.initial_sto_reg_iter) {
        sto_.setRegularization(solver_options_.initial_sto_reg);
      }
      else {
        sto_.setRegularization(0);
      }
      solver_statistics_.ts.emplace_back(contact_sequence_->eventTimes());
    } 
    updateSolution(t, q, v);
    const double kkt_error = KKTError();
    solver_statistics_.kkt_error.push_back(kkt_error); 
    if ((ocp_.sto_cost && ocp_.sto_constraints) && (kkt_error < solver_options_.kkt_tol_mesh)) {
      if (time_discretization_.dt_max() > solver_options_.max_dt_mesh) {
        if (solver_options_.enable_solution_interpolation) {
          solution_interpolator_.store(time_discretization_, s_);
        }
        discretize(t);
        if (solver_options_.enable_solution_interpolation) {
          solution_interpolator_.interpolateEventBased(robots_[0], time_discretization_, s_);
        }
        dms_.initConstraints(robots_, time_discretization_, s_);
        inner_iter = 0;
        solver_statistics_.mesh_refinement_iter.push_back(iter+1); 
      }
      else if (kkt_error < solver_options_.kkt_tol) {
        solver_statistics_.convergence = true;
        solver_statistics_.iter = iter+1;
        break;
      }
    }
    else if (kkt_error < solver_options_.kkt_tol) {
      solver_statistics_.convergence = true;
      solver_statistics_.iter = iter+1;
      break;
    }
  }
  if (!solver_statistics_.convergence) {
    solver_statistics_.iter = solver_options_.max_iter;
  }
  if (solver_options_.enable_solution_interpolation) {
    solution_interpolator_.store(time_discretization_, s_);
  }
  if (solver_options_.enable_benchmark) {
    timer_.tock();
    solver_statistics_.cpu_time = timer_.ms();
  }
}


const SolverStatistics& OCPSolver::getSolverStatistics() const {
  return solver_statistics_;
}


const Solution& OCPSolver::getSolution() const {
  return s_;
}


const SplitSolution& OCPSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage < time_discretization_.size());
  return s_[stage];
}


std::vector<Eigen::VectorXd> OCPSolver::getSolution(
    const std::string& name, const std::string& option) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      sol.push_back(s_[i].q);
    }
  }
  else if (name == "v") {
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      sol.push_back(s_[i].v);
    }
  }
  else if (name == "a") {
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      sol.push_back(s_[i].a);
    }
  }
  else if (name == "f" && option == "WORLD") {
    Robot robot = robots_[0];
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
      robot.updateFrameKinematics(s_[i].q);
      for (int j=0; j<robot.maxNumContacts(); ++j) {
        if (s_[i].isContactActive(j)) {
          const int contact_frame = robot.contactFrames()[j];
          robot.transformFromLocalToWorld(contact_frame, s_[i].f[j].template head<3>(),
                                          f.template segment<3>(3*j));
        }
      }
      sol.push_back(f);
    }
  }
  else if (name == "f") {
    Robot robot = robots_[0];
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
      for (int j=0; j<robot.maxNumContacts(); ++j) {
        if (s_[i].isContactActive(j)) {
          f.template segment<3>(3*j) = s_[i].f[j].template head<3>();
        }
      }
      sol.push_back(f);
    }
  }
  else if (name == "u") {
    for (int i=0; i<time_discretization_.size(); ++i) {
      if (time_discretization_[i].type == GridType::Impulse) continue;
      sol.push_back(s_[i].u);
    }
  }
  return sol;
}


const aligned_vector<LQRPolicy>& OCPSolver::getLQRPolicy() const {
  return riccati_recursion_.getLQRPolicy();
}


const RiccatiFactorization& OCPSolver::getRiccatiFactorization() const {
  return riccati_factorization_;
}


void OCPSolver::setSolution(const Solution& s) {
  s_ = s;
}


void OCPSolver::setSolution(const std::string& name, 
                            const Eigen::VectorXd& value) {
  if (name == "q") {
    if (value.size() != robots_[0].dimq()) {
      throw std::out_of_range(
          "[OCPSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
    }
    for (auto& e : s_) { e.q = value; }
  }
  else if (name == "v") {
    if (value.size() != robots_[0].dimv()) {
      throw std::out_of_range(
          "[OCPSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
    }
    for (auto& e : s_) { e.v = value; }
  }
  else if (name == "a") {
    if (value.size() != robots_[0].dimv()) {
      throw std::out_of_range(
          "[OCPSolver] invalid argument: a.size() must be " + std::to_string(robots_[0].dimv()) + "!");
    }
    for (auto& e : s_) { e.a = value; }
  }
  else if (name == "f") {
    if (value.size() == 6) {
      for (auto& e : s_) { 
        for (auto& ef : e.f) { ef = value.template head<6>(); } 
        e.set_f_stack(); 
      }
    }
    else if (value.size() == 3) {
      for (auto& e : s_) { 
        for (auto& ef : e.f) { ef.template head<3>() = value.template head<3>(); } 
        e.set_f_stack(); 
      }
    }
    else {
      throw std::out_of_range("[OCPSolver] invalid argument: f.size() must be 3 or 6!");
    }
  }
  else if (name == "u") {
    if (value.size() != robots_[0].dimu()) {
      throw std::out_of_range(
          "[OCPSolver] invalid argument: u.size() must be " + std::to_string(robots_[0].dimu()) + "!");
    }
    for (auto& e : s_) { e.u = value; }
  }
  else {
    throw std::invalid_argument("[OCPSolver] invalid arugment: name must be q, v, a, f, or u!");
  }
}


double OCPSolver::KKTError(const double t, const Eigen::VectorXd& q, 
                           const Eigen::VectorXd& v) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[OCPSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[OCPSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
  }
  resizeData();
  dms_.evalKKT(robots_, time_discretization_, q, v, s_, kkt_matrix_, kkt_residual_);
  sto_.evalKKT(time_discretization_, kkt_matrix_, kkt_residual_);
  return KKTError();
}


double OCPSolver::KKTError() const {
  return std::sqrt(dms_.getEval().kkt_error + sto_.getEval().kkt_error);
}


const TimeDiscretization& OCPSolver::getTimeDiscretization() const {
  return time_discretization_;
}


void OCPSolver::setRobotProperties(const RobotProperties& properties) {
  for (auto& e : robots_) {
    e.setRobotProperties(properties);
  }
}


template <typename T>
void conservativeReserve(const TimeDiscretization& time_discretization, 
                         aligned_vector<T>& data) {
  while (data.size() < time_discretization.size()) {
    data.push_back(data.back());
  }
}


void OCPSolver::resizeData() {
  conservativeReserve(time_discretization_, kkt_matrix_);
  conservativeReserve(time_discretization_, kkt_residual_);
  conservativeReserve(time_discretization_, s_);
  conservativeReserve(time_discretization_, d_);
  conservativeReserve(time_discretization_, riccati_factorization_);
  for (int i=0; i<time_discretization_.size(); ++i) {
    const auto& grid = time_discretization_[i];
    if (grid.type == GridType::Intermediate || grid.type == GridType::Lift) {
      s_[i].setContactStatus(contact_sequence_->contactStatus(grid.contact_phase));
      s_[i].set_f_stack();
    }
    else if (grid.type == GridType::Impulse) {
      s_[i].setContactStatus(contact_sequence_->impulseStatus(grid.impulse_index));
      s_[i].set_f_stack();
    }
    if (grid.switching_constraint) {
      const auto& grid_next_next = time_discretization_.grid(i+2);
      s_[i].setSwitchingConstraintDimension(contact_sequence_->impulseStatus(grid_next_next.impulse_index).dimf());
    }
    else {
      s_[i].setSwitchingConstraintDimension(0);
    }
  }
  dms_.resizeData(time_discretization_);
  riccati_recursion_.resizeData(time_discretization_);
  line_search_.resizeData(time_discretization_);
}


void OCPSolver::disp(std::ostream& os) const {
  os << ocp_ << std::endl;
}


std::ostream& operator<<(std::ostream& os, const OCPSolver& ocp_solver) {
  ocp_solver.disp(os);
  return os;
}

} // namespace robotoc