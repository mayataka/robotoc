#include "idocp/unocp/unparnmpc_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>
#include <fstream>


namespace idocp {

UnParNMPCSolver::UnParNMPCSolver(const Robot& robot, 
                                 const std::shared_ptr<CostFunction>& cost, 
                                 const std::shared_ptr<Constraints>& constraints, 
                                 const double T, const int N, const int num_proc)
  : robots_(num_proc, robot),
    parnmpc_(robot, cost, constraints, N),
    backward_correction_(robot, T, N, num_proc),
    unkkt_matrix_(N, SplitUnKKTMatrix(robot)),
    unkkt_residual_(N, SplitUnKKTResidual(robot)),
    s_(N, SplitSolution(robot)),
    d_(N, SplitDirection(robot)),
    N_(N),
    num_proc_(num_proc),
    T_(T),
    dtau_(T/N),
    kkt_error_(Eigen::VectorXd::Zero(N)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


UnParNMPCSolver::UnParNMPCSolver() {
}


UnParNMPCSolver::~UnParNMPCSolver() {
}


void UnParNMPCSolver::initConstraints() {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    parnmpc_[i].initConstraints(robots_[omp_get_thread_num()], i+1, s_[i]);
  }
}


void UnParNMPCSolver::initBackwardCorrection(const double t) {
  backward_correction_.initAuxMat(robots_, parnmpc_, t, s_);
}


void UnParNMPCSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  backward_correction_.coarseUpdate(robots_, parnmpc_, t, q, v, unkkt_matrix_,
                                    unkkt_residual_, s_, d_);
  backward_correction_.backwardCorrection(robots_, parnmpc_, s_, d_);
  const double primal_step_size = backward_correction_.primalStepSize();
  const double dual_step_size   = backward_correction_.dualStepSize();
  if (use_line_search) {
    // TODO: add filter line search method to choose primal_step_size
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      parnmpc_[i].updatePrimal(robots_[omp_get_thread_num()], primal_step_size, 
                               d_[i], s_[i]);
      parnmpc_[i].updateDual(dual_step_size);
    }
    else {
      parnmpc_.terminal.updatePrimal(robots_[omp_get_thread_num()],  
                                     primal_step_size, d_[i], s_[i]);
      parnmpc_.terminal.updateDual(dual_step_size);
    }
  }
} 


const SplitSolution& UnParNMPCSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


void UnParNMPCSolver::getStateFeedbackGain(const int time_stage, 
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


bool UnParNMPCSolver::setStateTrajectory(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  for (auto& e : s_) {
    e.v = v;
    e.q = q;
  }
  initConstraints();
  return isCurrentSolutionFeasible();
}


void UnParNMPCSolver::clearLineSearchFilter() {
  // filter_.clear();
}


double UnParNMPCSolver::KKTError() {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      kkt_error_.coeffRef(i) = parnmpc_[i].squaredNormKKTResidual(dtau_);
    }
    else {
      kkt_error_.coeffRef(N_-1) = parnmpc_.terminal.squaredNormKKTResidual(dtau_);
    }
  }
  return std::sqrt(kkt_error_.sum());
}


void UnParNMPCSolver::computeKKTResidual(const double t, 
                                         const Eigen::VectorXd& q, 
                                         const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc_[0].computeKKTResidual(robots_[omp_get_thread_num()], t+dtau_,   
                                     dtau_, q, v, s_[0], s_[1]);
    }
    else if (i < N_-1) {
      parnmpc_[i].computeKKTResidual(robots_[omp_get_thread_num()], t+(i+1)*dtau_,   
                                     dtau_, s_[i-1].q, s_[i-1].v, s_[i], s_[i+1]);
    }
    else {
      parnmpc_.terminal.computeKKTResidual(robots_[omp_get_thread_num()], t+T_, 
                                           dtau_, s_[i-1].q, s_[i-1].v, s_[i]);
    }
  }
}


bool UnParNMPCSolver::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = parnmpc_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time stage " << i << std::endl;
      return false;
    }
  }
  return true;
}


Robot UnParNMPCSolver::createRobot() const {
  return robots_[0];
}


std::vector<Eigen::VectorXd> UnParNMPCSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<N_; ++i) {
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


void UnParNMPCSolver::printSolution(const std::string& name, 
                                    const std::vector<int> frames) const {
  if (name == "all") {
    for (int i=0; i<N_; ++i) {
      std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
      std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
      std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
      std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    }
  }
  if (name == "q") {
    for (int i=0; i<N_; ++i) {
      std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    }
  }
  if (name == "v") {
    for (int i=0; i<N_; ++i) {
      std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    }
  }
  if (name == "a") {
    for (int i=0; i<N_; ++i) {
      std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    }
  }
  if (name == "u") {
    for (int i=0; i<N_; ++i) {
      std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    }
  }
  if (name == "end-effector") {
    Robot robot = robots_[0];
    for (int i=0; i<N_; ++i) {
      robot.updateFrameKinematics(s_[i].q);
      for (const auto e : frames) {
      std::cout << "end-effector[" << i << "][" << e << "] = " 
                << robot.framePosition(e).transpose() << std::endl;
      }
    }
  }
}


void UnParNMPCSolver::saveSolution(const std::string& path_to_file, 
                                   const std::string& name) const {
  std::ofstream file(path_to_file);
  if (name == "q") {
    const int dimq = robots_[0].dimq();
    for (int i=0; i<N_; ++i) {
      for (int j=0; j<dimq; ++j) {
        file << s_[i].q.coeff(j) << " ";
      }
      file << "\n";
    }
  }
  if (name == "v") {
    const int dimv = robots_[0].dimv();
    for (int i=0; i<N_; ++i) {
      for (int j=0; j<dimv; ++j) {
        file << s_[i].v.coeff(j) << " ";
      }
      file << "\n";
    }
  }
  if (name == "a") {
    const int dimv = robots_[0].dimv();
    for (int i=0; i<N_; ++i) {
      for (int j=0; j<dimv; ++j) {
        file << s_[i].a.coeff(j) << " ";
      }
      file << "\n";
    }
  }
  if (name == "u") {
    const int dimu = robots_[0].dimu();
    for (int i=0; i<N_; ++i) {
      for (int j=0; j<dimu; ++j) {
        file << s_[i].u.coeff(j) << " ";
      }
      file << "\n";
    }
  }
  file.close();
}

} // namespace idocp