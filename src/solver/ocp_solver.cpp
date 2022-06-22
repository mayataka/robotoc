#include "robotoc/solver/ocp_solver.hpp"

#include <stdexcept>
#include <cassert>
#include <algorithm>


namespace robotoc {

OCPSolver::OCPSolver(const OCP& ocp, 
                     const SolverOptions& solver_options, const int nthreads)
  : robots_(nthreads, ocp.robot()),
    contact_sequence_(ocp.contact_sequence()),
    dms_(nthreads),
    sto_(ocp),
    riccati_recursion_(ocp, nthreads, solver_options.max_dts_riccati),
    line_search_(ocp, nthreads),
    ocp_(ocp),
    riccati_factorization_(ocp.robot(), ocp.N(), ocp.maxNumEachDiscreteEvents()),
    kkt_matrix_(ocp.robot(), ocp.N(), ocp.maxNumEachDiscreteEvents()),
    kkt_residual_(ocp.robot(), ocp.N(), ocp.maxNumEachDiscreteEvents()),
    s_(ocp.robot(), ocp.N(), ocp.maxNumEachDiscreteEvents()),
    d_(ocp.robot(), ocp.N(), ocp.maxNumEachDiscreteEvents()),
    solver_options_(solver_options),
    solver_statistics_() {
  try {
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  for (auto& e : s_.data)    { ocp.robot().normalizeConfiguration(e.q); }
  for (auto& e : s_.impulse) { ocp.robot().normalizeConfiguration(e.q); }
  for (auto& e : s_.aux)     { ocp.robot().normalizeConfiguration(e.q); }
  for (auto& e : s_.lift)    { ocp.robot().normalizeConfiguration(e.q); }
}


OCPSolver::OCPSolver() {
}


OCPSolver::~OCPSolver() {
}


void OCPSolver::setSolverOptions(const SolverOptions& solver_options) {
  solver_options_ = solver_options;
  riccati_recursion_.setRegularization(solver_options.max_dts_riccati);
}


void OCPSolver::meshRefinement(const double t) {
  ocp_.meshRefinement(t);
  if (ocp_.discrete().discretizationMethod() == DiscretizationMethod::PhaseBased) {
    discretizeSolution();
    dms_.initConstraints(ocp_, robots_, contact_sequence_, s_);
    sto_.initConstraints(ocp_);
  }
}


void OCPSolver::initConstraints(const double t) {
  ocp_.discretize(t);
  discretizeSolution();
  dms_.initConstraints(ocp_, robots_, contact_sequence_, s_);
  sto_.initConstraints(ocp_);
}


void OCPSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  ocp_.discretize(t);
  discretizeSolution();
  dms_.computeKKTSystem(ocp_, robots_, contact_sequence_, q, v, s_, 
                        kkt_matrix_, kkt_residual_);
  sto_.computeKKTSystem(ocp_, kkt_matrix_, kkt_residual_);
  sto_.applyRegularization(ocp_, kkt_matrix_);
  riccati_recursion_.backwardRiccatiRecursion(ocp_, kkt_matrix_, kkt_residual_, 
                                              riccati_factorization_);
  dms_.computeInitialStateDirection(ocp_, robots_, q, v, s_, d_);
  riccati_recursion_.forwardRiccatiRecursion(ocp_, kkt_matrix_, kkt_residual_, d_);
  riccati_recursion_.computeDirection(ocp_, contact_sequence_, 
                                      riccati_factorization_, d_);
  sto_.computeDirection(ocp_, d_);
  double primal_step_size = std::min(riccati_recursion_.maxPrimalStepSize(), 
                                     sto_.maxPrimalStepSize());
  const double dual_step_size = std::min(riccati_recursion_.maxDualStepSize(),
                                         sto_.maxDualStepSize());
  if (solver_options_.enable_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(ocp_, robots_, 
                                                    contact_sequence_, 
                                                    q, v, s_, d_, 
                                                    max_primal_step_size);
  }
  solver_statistics_.primal_step_size.push_back(primal_step_size);
  solver_statistics_.dual_step_size.push_back(dual_step_size);
  dms_.integrateSolution(ocp_, robots_, primal_step_size, dual_step_size, 
                         kkt_matrix_, d_, s_);
  sto_.integrateSolution(ocp_, contact_sequence_, primal_step_size, 
                         dual_step_size, d_);
} 


void OCPSolver::solve(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, const bool init_solver) {
  if (solver_options_.enable_benchmark) {
    timer_.tick();
  }
  if (init_solver) {
    meshRefinement(t);
    initConstraints(t);
    line_search_.clearFilter();
  }
  solver_statistics_.clear(); 
  int inner_iter = 0;
  for (int iter=0; iter<solver_options_.max_iter; ++iter, ++inner_iter) {
    if (ocp_.isSTOEnabled()) {
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
    if (ocp_.isSTOEnabled() && (kkt_error < solver_options_.kkt_tol_mesh)) {
      if (ocp_.discrete().dt_max() > solver_options_.max_dt_mesh) {
        meshRefinement(t);
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
  assert(stage <= ocp_.discrete().N());
  return s_[stage];
}


std::vector<Eigen::VectorXd> OCPSolver::getSolution(
    const std::string& name, const std::string& option) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<=ocp_.discrete().N(); ++i) {
      sol.push_back(s_[i].q);
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        sol.push_back(s_.aux[impulse_index].q);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        sol.push_back(s_.lift[lift_index].q);
      }
    }
  }
  else if (name == "v") {
    for (int i=0; i<=ocp_.discrete().N(); ++i) {
      sol.push_back(s_[i].v);
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        sol.push_back(s_.aux[impulse_index].v);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        sol.push_back(s_.lift[lift_index].v);
      }
    }
  }
  else if (name == "a") {
    for (int i=0; i<ocp_.discrete().N(); ++i) {
      sol.push_back(s_[i].a);
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        sol.push_back(s_.aux[impulse_index].a);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        sol.push_back(s_.lift[lift_index].a);
      }
    }
  }
  else if (name == "f" && option == "WORLD") {
    Robot robot = robots_[0];
    for (int i=0; i<ocp_.discrete().N(); ++i) {
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
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
        robot.updateFrameKinematics(s_.aux[impulse_index].q);
        for (int j=0; j<robot.maxNumContacts(); ++j) {
          if (s_.aux[impulse_index].isContactActive(j)) {
            const int contact_frame = robot.contactFrames()[j];
            robot.transformFromLocalToWorld(contact_frame, 
                                            s_.aux[impulse_index].f[j].template head<3>(),
                                            f.template segment<3>(3*j));
          }
        }
        sol.push_back(f);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
        robot.updateFrameKinematics(s_.lift[lift_index].q);
        for (int j=0; j<robot.maxNumContacts(); ++j) {
          if (s_.lift[lift_index].isContactActive(j)) {
            const int contact_frame = robot.contactFrames()[j];
            robot.transformFromLocalToWorld(contact_frame, 
                                            s_.lift[lift_index].f[j].template head<3>(), 
                                            f.template segment<3>(3*j));
          }
        }
        sol.push_back(f);
      }
    }
  }
  else if (name == "f") {
    Robot robot = robots_[0];
    for (int i=0; i<ocp_.discrete().N(); ++i) {
      Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
      for (int j=0; j<robot.maxNumContacts(); ++j) {
        if (s_[i].isContactActive(j)) {
          f.template segment<3>(3*j) = s_[i].f[j].template head<3>();
        }
      }
      sol.push_back(f);
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
        for (int j=0; j<robot.maxNumContacts(); ++j) {
          if (s_.aux[impulse_index].isContactActive(j)) {
            f.template segment<3>(3*j) = s_.aux[impulse_index].f[j].template head<3>();
          }
        }
        sol.push_back(f);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        Eigen::VectorXd f(Eigen::VectorXd::Zero(robot.max_dimf()));
        for (int j=0; j<robot.maxNumContacts(); ++j) {
          if (s_.lift[lift_index].isContactActive(j)) {
            f.template segment<3>(3*j) = s_.lift[lift_index].f[j].template head<3>();
          }
        }
        sol.push_back(f);
      }
    }
  }
  else if (name == "u") {
    for (int i=0; i<ocp_.discrete().N(); ++i) {
      sol.push_back(s_[i].u);
      if (ocp_.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_.discrete().impulseIndexAfterTimeStage(i);
        sol.push_back(s_.aux[impulse_index].u);
      }
      else if (ocp_.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_.discrete().liftIndexAfterTimeStage(i);
        sol.push_back(s_.lift[lift_index].u);
      }
    }
  }
  else if (name == "ts") {
    const int num_events = ocp_.discrete().N_impulse()+ocp_.discrete().N_lift();
    int impulse_index = 0;
    int lift_index = 0;
    Eigen::VectorXd ts(num_events);
    for (int event_index=0; event_index<num_events; ++event_index) {
      if (ocp_.discrete().eventType(event_index) == DiscreteEventType::Impulse) {
        ts.coeffRef(event_index) = contact_sequence_->impulseTime(impulse_index);
        ++impulse_index;
      }
      else {
        ts.coeffRef(event_index) = contact_sequence_->liftTime(lift_index);
        ++lift_index;
      }
    }
    sol.push_back(ts);
  }
  return sol;
}


const hybrid_container<LQRPolicy>& OCPSolver::getLQRPolicy() const {
  return riccati_recursion_.getLQRPolicy();
}


const RiccatiFactorization& OCPSolver::getRiccatiFactorization() const {
  return riccati_factorization_;
}


void OCPSolver::setSolution(const Solution& s) {
  assert(s.data.size() == s_.data.size());
  assert(s.lift.size() == s_.lift.size());
  assert(s.aux.size() == s_.aux.size());
  assert(s.impulse.size() == s_.impulse.size());
  s_ = s;
}


void OCPSolver::setSolution(const std::string& name, 
                            const Eigen::VectorXd& value) {
  try {
    if (name == "q") {
      if (value.size() != robots_[0].dimq()) {
        throw std::out_of_range(
            "invalid value: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
      }
      for (auto& e : s_.data)    { e.q = value; }
      for (auto& e : s_.impulse) { e.q = value; }
      for (auto& e : s_.aux)     { e.q = value; }
      for (auto& e : s_.lift)    { e.q = value; }
    }
    else if (name == "v") {
      if (value.size() != robots_[0].dimv()) {
        throw std::out_of_range(
            "invalid value: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
      }
      for (auto& e : s_.data)    { e.v = value; }
      for (auto& e : s_.impulse) { e.v = value; }
      for (auto& e : s_.aux)     { e.v = value; }
      for (auto& e : s_.lift)    { e.v = value; }
    }
    else if (name == "a") {
      if (value.size() != robots_[0].dimv()) {
        throw std::out_of_range(
            "invalid value: a.size() must be " + std::to_string(robots_[0].dimv()) + "!");
      }
      for (auto& e : s_.data)    { e.a  = value; }
      for (auto& e : s_.impulse) { e.dv = value; }
      for (auto& e : s_.aux)     { e.a  = value; }
      for (auto& e : s_.lift)    { e.a  = value; }
    }
    else if (name == "f") {
      if (value.size() == 6) {
        for (auto& e : s_.data) { 
          for (auto& ef : e.f) { ef = value.template head<6>(); } 
          e.set_f_stack(); 
        }
        for (auto& e : s_.aux) { 
          for (auto& ef : e.f) { ef = value.template head<6>(); } 
          e.set_f_stack(); 
        }
        for (auto& e : s_.lift) { 
          for (auto& ef : e.f) { ef = value.template head<6>(); } 
          e.set_f_stack(); 
        }
      }
      else if (value.size() == 3) {
        for (auto& e : s_.data) { 
          for (auto& ef : e.f) { ef.template head<3>() = value.template head<3>(); } 
          e.set_f_stack(); 
        }
        for (auto& e : s_.aux) { 
          for (auto& ef : e.f) { ef.template head<3>() = value.template head<3>(); } 
          e.set_f_stack(); 
        }
        for (auto& e : s_.lift) { 
          for (auto& ef : e.f) { ef.template head<3>() = value.template head<3>(); } 
          e.set_f_stack(); 
        }
      }
      else {
        throw std::out_of_range("invalid value: f.size() must be 3 or 6!");
      }
    }
    else if (name == "lmd") {
      if (value.size() != 6) {
        for (auto& e : s_.impulse) { 
          for (auto& ef : e.f) { ef = value.template head<6>(); } 
          e.set_f_stack(); 
        }
      }
      else if (value.size() != 3) {
        for (auto& e : s_.impulse) { 
          for (auto& ef : e.f) { ef.template head<3>() = value.template head<3>(); } 
          e.set_f_stack(); 
        }
      }
      else {
        throw std::out_of_range("invalid value: f.size() must be 3 or 6!");
      }
    }
    else if (name == "u") {
      if (value.size() != robots_[0].dimu()) {
        throw std::out_of_range(
            "invalid value: u.size() must be " + std::to_string(robots_[0].dimu()) + "!");
      }
      for (auto& e : s_.data)    { e.u = value; }
      for (auto& e : s_.aux)     { e.u = value; }
      for (auto& e : s_.lift)    { e.u = value; }
    }
    else {
      throw std::invalid_argument("invalid arugment: name must be q, v, a, f, or u!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


void OCPSolver::extrapolateSolutionLastPhase(const double t) {
  const int num_discrete_events = contact_sequence_->numDiscreteEvents();
  if (num_discrete_events > 0) {
    ocp_.discretize(t);
    int time_stage_after_last_event;
    if (contact_sequence_->eventType(num_discrete_events-1) 
          == DiscreteEventType::Impulse) {
      time_stage_after_last_event 
          = ocp_.discrete().timeStageAfterImpulse(ocp_.discrete().N_impulse()-1);
    }
    else {
      time_stage_after_last_event 
          = ocp_.discrete().timeStageAfterLift(ocp_.discrete().N_lift()-1);
    }
    for (int i=time_stage_after_last_event; i<=ocp_.discrete().N(); ++i) {
      s_[i].copyPrimal(s_[time_stage_after_last_event-1]);
      s_[i].copyDual(s_[time_stage_after_last_event-1]);
      ocp_[i].initConstraints(ocp_[time_stage_after_last_event-1]);
    }
  }
}


void OCPSolver::extrapolateSolutionInitialPhase(const double t) {
  const int num_discrete_events = contact_sequence_->numDiscreteEvents();
  if (num_discrete_events > 0) {
    ocp_.discretize(t);
    int time_stage_before_initial_event;
    if (contact_sequence_->eventType(0) == DiscreteEventType::Impulse) {
      time_stage_before_initial_event 
          = ocp_.discrete().timeStageBeforeImpulse(0);
    }
    else {
      time_stage_before_initial_event 
          = ocp_.discrete().timeStageBeforeLift(0);
    }
    for (int i=0; i<=time_stage_before_initial_event; ++i) {
      s_[i].copyPrimal(s_[time_stage_before_initial_event+1]);
      s_[i].copyDual(s_[time_stage_before_initial_event+1]);
      ocp_[i].initConstraints(ocp_[time_stage_before_initial_event+1]);
    }
  }
}


double OCPSolver::KKTError(const double t, const Eigen::VectorXd& q, 
                           const Eigen::VectorXd& v) {
  ocp_.discretize(t);
  discretizeSolution();
  dms_.computeKKTResidual(ocp_, robots_, contact_sequence_, q, v, s_, 
                          kkt_matrix_, kkt_residual_);
  sto_.computeKKTResidual(ocp_, kkt_residual_);
  return KKTError();
}


double OCPSolver::KKTError() const {
  return std::sqrt(dms_.KKTError(ocp_, kkt_residual_) + sto_.KKTError());
}


double OCPSolver::cost(const bool include_cost_barrier) const {
  return dms_.totalCost(ocp_, include_cost_barrier);
}


bool OCPSolver::isCurrentSolutionFeasible(const bool verbose) {
  // ocp_.discretize(t);
  // discretizeSolution();
  return dms_.isFeasible(ocp_, robots_, contact_sequence_, s_);
}


const TimeDiscretization& OCPSolver::getTimeDiscretization() const {
  return ocp_.discrete();
}


void OCPSolver::discretizeSolution() {
  for (int i=0; i<=ocp_.discrete().N(); ++i) {
    s_[i].setContactStatus(
        contact_sequence_->contactStatus(ocp_.discrete().contactPhase(i)));
    s_[i].set_f_stack();
    s_[i].setImpulseStatus();
  }
  for (int i=0; i<ocp_.discrete().N_lift(); ++i) {
    s_.lift[i].setContactStatus(
        contact_sequence_->contactStatus(
            ocp_.discrete().contactPhaseAfterLift(i)));
    s_.lift[i].set_f_stack();
    s_.lift[i].setImpulseStatus();
  }
  for (int i=0; i<ocp_.discrete().N_impulse(); ++i) {
    s_.impulse[i].setImpulseStatus(contact_sequence_->impulseStatus(i));
    s_.impulse[i].set_f_stack();
    s_.aux[i].setContactStatus(
        contact_sequence_->contactStatus(
            ocp_.discrete().contactPhaseAfterImpulse(i)));
    s_.aux[i].set_f_stack();
    const int time_stage_before_impulse 
        = ocp_.discrete().timeStageBeforeImpulse(i);
    if (time_stage_before_impulse-1 >= 0) {
      s_[time_stage_before_impulse-1].setImpulseStatus(
          contact_sequence_->impulseStatus(i));
    }
  }
}


void OCPSolver::disp(std::ostream& os) const {
  os << ocp_ << std::endl;
}


std::ostream& operator<<(std::ostream& os, const OCPSolver& ocp_solver) {
  ocp_solver.disp(os);
  return os;
}

} // namespace robotoc