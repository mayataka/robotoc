#include "idocp/solver/unocp_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

UnOCPSolver::UnOCPSolver(const Robot& robot, 
                         const std::shared_ptr<CostFunction>& cost, 
                         const std::shared_ptr<Constraints>& constraints, 
                         const double T, const int N, const int nthreads)
  : robots_(nthreads, robot),
    ocp_(robot, cost, constraints, N),
    riccati_recursion_(robot, T, N),
    line_search_(robot, T, N, nthreads),
    terminal_kkt_matrix_(robot),
    terminal_kkt_residual_(robot),
    unkkt_matrix_(N, SplitUnKKTMatrix(robot)),
    unkkt_residual_(N, SplitUnKKTResidual(robot)),
    s_(N+1, SplitSolution(robot)),
    d_(N+1, SplitDirection(robot)),
    riccati_factorization_(N+1, SplitRiccatiFactorization(robot)),
    N_(N),
    nthreads_(nthreads),
    T_(T),
    dt_(T/N),
    primal_step_size_(Eigen::VectorXd::Zero(N)), 
    dual_step_size_(Eigen::VectorXd::Zero(N)), 
    kkt_error_(Eigen::VectorXd::Zero(N+1)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


UnOCPSolver::UnOCPSolver() {
}


UnOCPSolver::~UnOCPSolver() {
}


void UnOCPSolver::initConstraints() {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      ocp_[i].initConstraints(robots_[omp_get_thread_num()], i, s_[i]);
    }
    else {
      ocp_.terminal.initConstraints(robots_[omp_get_thread_num()], N_, s_[N_]);
    }
  }
}


void UnOCPSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const bool line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].linearizeOCP(robots_[omp_get_thread_num()], t, dt_, q, 
                           s_[0], s_[1], unkkt_matrix_[0], unkkt_residual_[0]);
    }
    else if (i < N_) {
      ocp_[i].linearizeOCP(robots_[omp_get_thread_num()], t+i*dt_, dt_, 
                           s_[i-1].q, s_[i], s_[i+1], 
                           unkkt_matrix_[i], unkkt_residual_[i]);
    }
    else {
      ocp_.terminal.linearizeOCP(robots_[omp_get_thread_num()], t+T_, 
                                 s_[N_-1].q, s_[N_], 
                                 terminal_kkt_matrix_, terminal_kkt_residual_);
    }
  }
  riccati_recursion_.backwardRiccatiRecursionTerminal(terminal_kkt_matrix_, 
                                                      terminal_kkt_residual_,
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(unkkt_matrix_, unkkt_residual_,
                                              riccati_factorization_);
  d_[0].dq() = q - s_[0].q;
  d_[0].dv() = v - s_[0].v;
  riccati_recursion_.forwardRiccatiRecursion(unkkt_residual_, d_);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    SplitUnRiccatiFactorizer::computeCostateDirection(riccati_factorization_[i], 
                                                      d_[i]);
    if (i < N_) {
      ocp_[i].computeCondensedDirection(robots_[omp_get_thread_num()], dt_,
                                        s_[i], d_[i]);
      primal_step_size_.coeffRef(i) = ocp_[i].maxPrimalStepSize();
      dual_step_size_.coeffRef(i)   = ocp_[i].maxDualStepSize();
    }
  }
  double primal_step_size = primal_step_size_.minCoeff();
  const double dual_step_size   = dual_step_size_.minCoeff();
  if (line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(ocp_, robots_, t, q, v, s_,
                                                    d_, max_primal_step_size);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      ocp_[i].updatePrimal(robots_[omp_get_thread_num()], primal_step_size, 
                           d_[i], s_[i]);
      ocp_[i].updateDual(dual_step_size);
    }
    else {
      ocp_.terminal.updatePrimal(robots_[omp_get_thread_num()],  
                                 primal_step_size, d_[N_], s_[N_]);
      ocp_.terminal.updateDual(dual_step_size);
    }
  }
} 


const SplitSolution& UnOCPSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


std::vector<Eigen::VectorXd> UnOCPSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<=N_; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<=N_; ++i) {
      sol.push_back(s_[i].v);
    }
  }
  if (name == "a") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].a);
    }
  }
  if (name == "u") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].u);
    }
  }
  return sol;
}


void UnOCPSolver::getStateFeedbackGain(const int time_stage, 
                                       Eigen::MatrixXd& Kq, 
                                       Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  // riccati_solver_.getStateFeedbackGain(time_stage, Kq, Kv);
}


void UnOCPSolver::setSolution(const std::string& name, 
                              const Eigen::VectorXd& value) {
  try {
    if (name == "q") {
      for (auto& e : s_) { e.q = value; }
    }
    else if (name == "v") {
      for (auto& e : s_) { e.v = value; }
    }
    else if (name == "a") {
      for (auto& e : s_) { e.a  = value; }
    }
    else if (name == "u") {
      for (auto& e : s_) { e.u = value; }
    }
    else {
      throw std::invalid_argument("invalid arugment: name must be q, v, a, or u!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


void UnOCPSolver::clearLineSearchFilter() {
  line_search_.clearFilter();
}


void UnOCPSolver::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].computeKKTResidual(robots_[omp_get_thread_num()], t, dt_,  
                                 q, s_[0], s_[1]);
    }
    else if (i < N_) {
      ocp_[i].computeKKTResidual(robots_[omp_get_thread_num()], t+i*dt_, dt_,  
                                 s_[i-1].q, s_[i], s_[i+1]);
    }
    else {
      ocp_.terminal.computeKKTResidual(robots_[omp_get_thread_num()], t+T_, 
                                       s_[N_-1].q, s_[N_], terminal_kkt_matrix_, 
                                       terminal_kkt_residual_);
    }
  }
}


double UnOCPSolver::KKTError() {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      kkt_error_.coeffRef(i) = ocp_[i].squaredNormKKTResidual(dt_);
    }
    else {
      kkt_error_.coeffRef(N_) 
          = ocp_.terminal.squaredNormKKTResidual(terminal_kkt_residual_);
    }
  }
  return std::sqrt(kkt_error_.sum());
}


bool UnOCPSolver::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = ocp_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time stage " << i << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace idocp