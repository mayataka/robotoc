#include "idocp/ocp/backward_correction.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

BackwardCorrection::BackwardCorrection(const Robot& robot, const int N, 
                                       const int max_num_impulse, 
                                       const int nthreads)
  : N_(N),
    max_num_impulse_(max_num_impulse),
    nthreads_(nthreads),
    N_all_(N),
    corrector_(robot, N, max_num_impulse),
    s_new_(robot, N, max_num_impulse),
    aux_mat_(N, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_impulse_(max_num_impulse, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_aux_(max_num_impulse, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_lift_(max_num_impulse, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N+3*max_num_impulse)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N+3*max_num_impulse)) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
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


BackwardCorrection::BackwardCorrection()
  : N_(0),
    max_num_impulse_(0),
    nthreads_(0),
    N_all_(0),
    corrector_(),
    s_new_(),
    aux_mat_(),
    aux_mat_impulse_(),
    aux_mat_aux_(),
    aux_mat_lift_(),
    primal_step_sizes_(),
    dual_step_sizes_() {
}


BackwardCorrection::~BackwardCorrection() {
}


void BackwardCorrection::initAuxMat(ParNMPC& parnmpc, 
                                    std::vector<Robot>& robots, const double t, 
                                    const Solution& s, KKTMatrix& kkt_matrix) {
  parnmpc.terminal.computeTerminalCostHessian(robots[0], 
                                              parnmpc.discrete().t(N_-1), 
                                              s[N_-1], kkt_matrix[N_-1]);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    aux_mat_[i] = kkt_matrix[N_-1].Qxx();
  }
}


void BackwardCorrection::coarseUpdate(ParNMPC& parnmpc, 
                                      std::vector<Robot>& robots, 
                                      const ContactSequence& contact_sequence, 
                                      const double t, const Eigen::VectorXd& q, 
                                      const Eigen::VectorXd& v,
                                      const Solution& s, KKTMatrix& kkt_matrix, 
                                      KKTResidual& kkt_residual) {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc[i].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
          parnmpc.discrete().t(i), parnmpc.discrete().dtau(i),
          q, v, s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                 parnmpc.discrete().dtau(i), aux_mat_[i+1], 
                                 kkt_matrix[i], kkt_residual[i], 
                                 s[i], s_new_[i]);
    }
    else if (i < N_-1) {
      parnmpc[i].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
          parnmpc.discrete().t(i), parnmpc.discrete().dtau(i),
          s[i-1].q, s[i-1].v, s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                 parnmpc.discrete().dtau(i), aux_mat_[i+1], 
                                 kkt_matrix[i], kkt_residual[i], 
                                 s[i], s_new_[i]);
    }
    else {
      parnmpc.terminal.linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)),  
          parnmpc.discrete().t(i), parnmpc.discrete().dtau(i),
          s[i-1].q, s[i-1].v, s[i], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                 parnmpc.discrete().dtau(i), kkt_matrix[i], 
                                 kkt_residual[i], s[i], s_new_[i]);
    }
  }
}


void BackwardCorrection::backwardCorrection(ParNMPC& parnmpc, 
                                            std::vector<Robot>& robots, 
                                            const KKTMatrix& kkt_matrix, 
                                            const KKTResidual& kkt_residual,
                                            const Solution& s, Direction& d) {
  for (int i=N_-2; i>=0; --i) {
    corrector_[i].backwardCorrectionSerial(s[i+1], s_new_[i+1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=N_-2; i>=0; --i) {
    corrector_[i].backwardCorrectionParallel(robots[omp_get_thread_num()], 
                                             s_new_[i]);
  }
  for (int i=1; i<N_; ++i) {
    corrector_[i].forwardCorrectionSerial(robots[0], s[i-1], s_new_[i-1], 
                                          s_new_[i]);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i > 0) {
      corrector_[i].forwardCorrectionParallel(s_new_[i]);
      aux_mat_[i] = - corrector_[i].auxMat();
    }
    SplitBackwardCorrection::computeDirection(robots[omp_get_thread_num()], 
                                              s[i], s_new_[i], d[i]);
    if (i < N_-1) {
      parnmpc[i].computeCondensedPrimalDirection(robots[omp_get_thread_num()], 
                                                 parnmpc.discrete().dtau(i), 
                                                 s[i], d[i]);
      parnmpc[i].computeCondensedDualDirection(robots[omp_get_thread_num()], 
                                               parnmpc.discrete().dtau(i), 
                                               kkt_matrix[i], kkt_residual[i], d[i]);
      primal_step_sizes_.coeffRef(i) = parnmpc[i].maxPrimalStepSize();
      dual_step_sizes_.coeffRef(i)   = parnmpc[i].maxDualStepSize();
    }
    else {
      parnmpc.terminal.computeCondensedPrimalDirection(
          robots[omp_get_thread_num()], parnmpc.discrete().dtau(i), s[i], d[i]);
      parnmpc.terminal.computeCondensedDualDirection(
          robots[omp_get_thread_num()], parnmpc.discrete().dtau(i), 
          kkt_matrix[i], kkt_residual[i], d[i]);
      primal_step_sizes_.coeffRef(i) = parnmpc.terminal.maxPrimalStepSize();
      dual_step_sizes_.coeffRef(i)   = parnmpc.terminal.maxDualStepSize();
    }
  }
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  N_all_ = N_ + 2 * N_impulse + N_lift;
}


double BackwardCorrection::primalStepSize() const {
  return primal_step_sizes_.head(N_all_).minCoeff();
}


double BackwardCorrection::dualStepSize() const {
  return dual_step_sizes_.head(N_all_).minCoeff();
}

} // namespace idocp