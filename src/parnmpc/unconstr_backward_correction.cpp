#include "idocp/parnmpc/unconstr_backward_correction.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

UnconstrBackwardCorrection::UnconstrBackwardCorrection(const Robot& robot, 
                                                       const double T, 
                                                       const int N, 
                                                       const int nthreads)
  : N_(N),
    nthreads_(nthreads),
    T_(T),
    dt_(T/N),
    corrector_(N, UnconstrSplitBackwardCorrection(robot)),
    s_new_(robot, N),
    aux_mat_(N, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)) {
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
}


UnconstrBackwardCorrection::UnconstrBackwardCorrection()
  : N_(0),
    nthreads_(0),
    T_(0),
    dt_(0),
    corrector_(),
    s_new_(),
    aux_mat_(),
    primal_step_sizes_(),
    dual_step_sizes_() {
}


UnconstrBackwardCorrection::~UnconstrBackwardCorrection() {
}


void UnconstrBackwardCorrection::initAuxMat(aligned_vector<Robot>& robots, 
                                            UnconstrParNMPC& parnmpc, 
                                            const double t, const Solution& s, 
                                            KKTMatrix& kkt_matrix,
                                            KKTResidual& kkt_residual) {
  parnmpc.terminal.computeTerminalCostHessian(robots[0], t+T_, s[N_-1], 
                                              kkt_matrix[0], kkt_residual[0]);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    aux_mat_[i] = kkt_matrix[0].Qxx;
  }
  kkt_matrix[0].setZero();
  kkt_residual[0].setZero();
}


void UnconstrBackwardCorrection::coarseUpdate(aligned_vector<Robot>& robots, 
                                              UnconstrParNMPC& parnmpc, 
                                              const double t, 
                                              const Eigen::VectorXd& q, 
                                              const Eigen::VectorXd& v,
                                              KKTMatrix& kkt_matrix, 
                                              KKTResidual& kkt_residual,
                                              const Solution& s) {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc[i].linearizeOCP(robots[omp_get_thread_num()], t+dt_, dt_, q, v, 
                              s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(aux_mat_[i+1], dt_, kkt_matrix[i], 
                                 kkt_residual[i], s[i], s_new_[i]);
    }
    else if (i < N_-1) {
      parnmpc[i].linearizeOCP(robots[omp_get_thread_num()], t+(i+1)*dt_, 
                              dt_, s[i-1].q, s[i-1].v, s[i], s[i+1], 
                              kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(aux_mat_[i+1], dt_, kkt_matrix[i], 
                                 kkt_residual[i], s[i], s_new_[i]);
    }
    else {
      parnmpc.terminal.linearizeOCP(robots[omp_get_thread_num()], t+T_, dt_, 
                                    s[i-1].q, s[i-1].v, s[i], 
                                    kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(dt_, kkt_matrix[i], kkt_residual[i], 
                                 s[i], s_new_[i]);
    }
  }
}


void UnconstrBackwardCorrection::backwardCorrection(UnconstrParNMPC& parnmpc, 
                                                    const Solution& s, 
                                                    const KKTMatrix& kkt_matrix, 
                                                    const KKTResidual& kkt_residual,
                                                    Direction& d) {
  for (int i=N_-2; i>=0; --i) {
    corrector_[i].backwardCorrectionSerial(s[i+1], s_new_[i+1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_) 
  for (int i=N_-2; i>=0; --i) {
    corrector_[i].backwardCorrectionParallel(s_new_[i]);
  }
  for (int i=1; i<N_; ++i) {
    corrector_[i].forwardCorrectionSerial(s[i-1], s_new_[i-1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i > 0) {
      corrector_[i].forwardCorrectionParallel(s_new_[i]);
      aux_mat_[i] = - corrector_[i].auxMat();
    }
    UnconstrSplitBackwardCorrection::computeDirection(s[i], s_new_[i], d[i]);
    if (i < N_-1) {
      parnmpc[i].expandPrimalAndDual(dt_, s[i], kkt_matrix[i], 
                                     kkt_residual[i], d[i]);
      primal_step_sizes_.coeffRef(i) = parnmpc[i].maxPrimalStepSize();
      dual_step_sizes_.coeffRef(i)   = parnmpc[i].maxDualStepSize();
    }
    else {
      parnmpc.terminal.expandPrimalAndDual(dt_, s[i], kkt_matrix[i], 
                                           kkt_residual[i], d[i]);
      primal_step_sizes_.coeffRef(i) = parnmpc.terminal.maxPrimalStepSize();
      dual_step_sizes_.coeffRef(i)   = parnmpc.terminal.maxDualStepSize();
    }
  }
}


double UnconstrBackwardCorrection::primalStepSize() const {
  return primal_step_sizes_.minCoeff();
}


double UnconstrBackwardCorrection::dualStepSize() const {
  return dual_step_sizes_.minCoeff();
}

} // namespace idocp