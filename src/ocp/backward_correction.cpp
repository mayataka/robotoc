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
    aux_mat_impulse_(max_num_impulse, 
                     Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_aux_(max_num_impulse, 
                 Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_lift_(max_num_impulse, 
                  Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
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
                                    std::vector<Robot>& robots, 
                                    const Solution& s, KKTMatrix& kkt_matrix) {
  parnmpc.terminal.computeTerminalCostHessian(robots[0], 
                                              parnmpc.discrete().t(N_-1), 
                                              s[N_-1], kkt_matrix[N_-1]);
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  N_all_ = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all_; ++i) {
    if (i < N_) {
      aux_mat_[i] = kkt_matrix[N_-1].Qxx();
    }
    else if (i < N_+N_impulse) {
      const int impulse_index = i - N_;
      aux_mat_impulse_[impulse_index] = kkt_matrix[N_-1].Qxx();
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index = i - N_ - N_impulse;
      aux_mat_aux_[impulse_index] = kkt_matrix[N_-1].Qxx();
    }
    else {
      const int lift_index = i - N_ - 2*N_impulse;
      aux_mat_lift_[lift_index] = kkt_matrix[N_-1].Qxx();
    }
  }
}


void BackwardCorrection::coarseUpdate(ParNMPC& parnmpc, 
                                      std::vector<Robot>& robots, 
                                      const ContactSequence& contact_sequence, 
                                      const Eigen::VectorXd& q, 
                                      const Eigen::VectorXd& v,
                                      const Solution& s, KKTMatrix& kkt_matrix, 
                                      KKTResidual& kkt_residual) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  N_all_ = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all_; ++i) {
    if (i < N_-1) {
      if (parnmpc.discrete().isTimeStageBeforeImpulse(i)) {
        const int impulse_index 
            = parnmpc.discrete().impulseIndexAfterTimeStage(i);
        parnmpc[i].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), s[i], s.aux[impulse_index], 
            kkt_matrix[i], kkt_residual[i]);
        corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                   parnmpc.discrete().dtau(i), 
                                   aux_mat_aux_[impulse_index], kkt_matrix[i], 
                                   kkt_residual[i], s[i], s_new_[i]);
      }
      else if (parnmpc.discrete().isTimeStageBeforeLift(i)) {
        const int lift_index = parnmpc.discrete().liftIndexAfterTimeStage(i);
        parnmpc[i].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), s[i], s.lift[lift_index], 
            kkt_matrix[i], kkt_residual[i]);
        corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                   parnmpc.discrete().dtau(i), 
                                   aux_mat_lift_[lift_index], kkt_matrix[i], 
                                   kkt_residual[i], s[i], s_new_[i]);
      }
      else {
        parnmpc[i].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), 
            s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
        corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                   parnmpc.discrete().dtau(i), 
                                   aux_mat_[i+1], kkt_matrix[i], 
                                   kkt_residual[i], s[i], s_new_[i]);
      }
    }
    else if (i == N_-1) {
      parnmpc.terminal.linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
          parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
          q_prev(parnmpc.discrete(), q, s, i), 
          v_prev(parnmpc.discrete(), v, s, i), 
          s[i], kkt_matrix[i], kkt_residual[i]);
      corrector_[i].coarseUpdate(robots[omp_get_thread_num()], 
                                 parnmpc.discrete().dtau(i), kkt_matrix[i], 
                                 kkt_residual[i], s[i], s_new_[i]);
    }
    else if (i < N_+N_impulse) {
      const int impulse_index = i - N_;
      const int time_stage_after_impulse 
          = parnmpc.discrete().timeStageAfterImpulse(impulse_index);
      parnmpc.impulse[impulse_index].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.impulseStatus(impulse_index), 
          parnmpc.discrete().t_impulse(impulse_index), s.aux[impulse_index].q, 
          s.aux[impulse_index].v, s.impulse[impulse_index], 
          s[time_stage_after_impulse], kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index]);
      corrector_.impulse[i].coarseUpdate(robots[omp_get_thread_num()], 
                                         aux_mat_[time_stage_after_impulse], 
                                         kkt_matrix.impulse[impulse_index], 
                                         kkt_residual.impulse[impulse_index], 
                                         s.impulse[impulse_index], 
                                         s_new_.impulse[impulse_index]);
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index  = i - (N_+N_impulse);
      const int time_stage_before_impulse 
          = parnmpc.discrete().timeStageBeforeImpulse(impulse_index);
      if (time_stage_before_impulse >= 0) {
        parnmpc.aux[impulse_index].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseBeforeImpulse(impulse_index)), 
            contact_sequence.impulseStatus(impulse_index),
            parnmpc.discrete().t_impulse(impulse_index), 
            parnmpc.discrete().dtau_aux(impulse_index), 
            s[time_stage_before_impulse].q, s[time_stage_before_impulse].v, 
            s.aux[impulse_index], s.impulse[impulse_index], 
            kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
      }
      else {
        assert(time_stage_before_impulse == -1);
        parnmpc.aux[impulse_index].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseBeforeImpulse(impulse_index)), 
            parnmpc.discrete().t_impulse(impulse_index), 
            parnmpc.discrete().dtau_aux(impulse_index), 
            q, v, s.aux[impulse_index], s.impulse[impulse_index], 
            kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
      }
      corrector_.aux[i].coarseUpdate(robots[omp_get_thread_num()], 
                                     parnmpc.discrete().dtau_aux(impulse_index), 
                                     aux_mat_impulse_[impulse_index], 
                                     kkt_matrix.aux[impulse_index], 
                                     kkt_residual.aux[impulse_index], 
                                     s.aux[impulse_index], 
                                     s_new_.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+2*N_impulse);
      const int time_stage_before_lift
          = parnmpc.discrete().timeStageBeforeLift(lift_index);
      if (time_stage_before_lift >= 0) {
        parnmpc.lift[lift_index].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseBeforeLift(lift_index)), 
            parnmpc.discrete().t_lift(lift_index), 
            parnmpc.discrete().dtau_lift(lift_index), 
            s[time_stage_before_lift].q, s[time_stage_before_lift].v, 
            s.lift[lift_index], s[time_stage_before_lift+1], 
            kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
      }
      else {
        assert(time_stage_before_lift == -1);
        parnmpc.lift[lift_index].linearizeOCP(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseBeforeLift(lift_index)), 
            parnmpc.discrete().t_lift(lift_index), 
            parnmpc.discrete().dtau_lift(lift_index), 
            q, v, s.lift[lift_index], s[time_stage_before_lift+1], 
            kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
      }
      corrector_.lift[i].coarseUpdate(robots[omp_get_thread_num()], 
                                      parnmpc.discrete().dtau_lift(lift_index), 
                                      aux_mat_[time_stage_before_lift+1], 
                                      kkt_matrix.lift[lift_index], 
                                      kkt_residual.lift[lift_index], 
                                      s.lift[lift_index], 
                                      s_new_.lift[lift_index]);
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
}


double BackwardCorrection::primalStepSize() const {
  return primal_step_sizes_.head(N_all_).minCoeff();
}


double BackwardCorrection::dualStepSize() const {
  return dual_step_sizes_.head(N_all_).minCoeff();
}

} // namespace idocp