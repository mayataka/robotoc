#include "idocp/ocp/parnmpc_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp{

ParNMPCLinearizer::ParNMPCLinearizer(const int N, const int max_num_impulse, 
                                     const int nthreads) 
  : N_(N),
    max_num_impulse_(max_num_impulse),
    nthreads_(nthreads),
    kkt_error_(Eigen::VectorXd::Zero(N+3*max_num_impulse)) {
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


ParNMPCLinearizer::ParNMPCLinearizer()
  : N_(0),
    max_num_impulse_(0),
    nthreads_(0),
    kkt_error_() {
}


ParNMPCLinearizer::~ParNMPCLinearizer() {
}


void ParNMPCLinearizer::initConstraints(ParNMPC& parnmpc, 
                                        std::vector<Robot>& robots, 
                                        const ContactSequence& contact_sequence, 
                                        const Solution& s) const {
  const int N_impulse = max_num_impulse_;
  const int N_lift = max_num_impulse_;
  const int N_all = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_-1) {
      parnmpc[i].initConstraints(robots[omp_get_thread_num()], i+1, s[i]);
    }
    else if (i == N_-1) {
      parnmpc.terminal.initConstraints(robots[omp_get_thread_num()], N_, s[N_]);
    }
    else if (i < N_+N_impulse) {
      const int impulse_index  = i - N_;
      parnmpc.impulse[impulse_index].initConstraints(
          robots[omp_get_thread_num()], s.impulse[impulse_index]);
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index  = i - (N_+N_impulse);
      parnmpc.aux[impulse_index].initConstraints(robots[omp_get_thread_num()], 
                                                 0, s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+2*N_impulse);
      parnmpc.lift[lift_index].initConstraints(robots[omp_get_thread_num()], 
                                               0, s.lift[lift_index]);
    }
  }
}


void ParNMPCLinearizer::computeKKTResidual(
    ParNMPC& parnmpc, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  const int N_all = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_-1) {
      if (parnmpc.discrete().isTimeStageBeforeImpulse(i)) {
        parnmpc[i].computeKKTResidual(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), 
            s[i], s.impulse[parnmpc.discrete().impulseIndexAfterTimeStage(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else if (parnmpc.discrete().isTimeStageBeforeLift(i)) {
        parnmpc[i].computeKKTResidual(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), 
            s[i], s.lift[parnmpc.discrete().liftIndexAfterTimeStage(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else {
        parnmpc[i].computeKKTResidual(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
            parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
            q_prev(parnmpc.discrete(), q, s, i), 
            v_prev(parnmpc.discrete(), v, s, i), 
            s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N_-1) {
      parnmpc.terminal.computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(parnmpc.discrete().contactPhase(i)), 
          parnmpc.discrete().t(i), parnmpc.discrete().dtau(i), 
          q_prev(parnmpc.discrete(), q, s, i), 
          v_prev(parnmpc.discrete(), v, s, i), 
          s[i], kkt_matrix[i], kkt_residual[i]);
    }
    else if (i < N_+N_impulse) {
      const int impulse_index  = i - N_;
      const int time_stage_after_impulse 
          = parnmpc.discrete().timeStageBeforeImpulse(impulse_index);
      parnmpc.impulse[impulse_index].computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.impulseStatus(impulse_index), 
          parnmpc.discrete().t_impulse(impulse_index), s.aux[impulse_index].q, 
          s.aux[impulse_index].v, s.impulse[impulse_index], 
          s[time_stage_after_impulse], kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index  = i - (N_+N_impulse);
      const int time_stage_before_impulse 
          = parnmpc.discrete().timeStageBeforeImpulse(impulse_index);
      if (time_stage_before_impulse >= 0) {
        parnmpc.aux[impulse_index].computeKKTResidual(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseAfterImpulse(impulse_index)), 
            contact_sequence.impulseStatus(impulse_index),
            parnmpc.discrete().t_impulse(impulse_index), 
            parnmpc.discrete().dtau_aux(impulse_index), 
            s[time_stage_before_impulse].q, s[time_stage_before_impulse].v, 
            s.aux[impulse_index], s.impulse[impulse_index], 
            kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
      }
      else {
        parnmpc.aux[impulse_index].computeKKTResidual(
            robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(
                parnmpc.discrete().contactPhaseAfterImpulse(impulse_index)), 
            parnmpc.discrete().t_impulse(impulse_index), 
            parnmpc.discrete().dtau_aux(impulse_index), 
            s[time_stage_before_impulse].q, s[time_stage_before_impulse].v, 
            s.aux[impulse_index], s.impulse[impulse_index], 
            kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
      }
    }
    else {
      const int lift_index = i - (N_+2*N_impulse);
      const int time_stage_before_lift
          = parnmpc.discrete().timeStageBeforeLift(lift_index);
      parnmpc.lift[lift_index].computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              parnmpc.discrete().contactPhaseAfterLift(lift_index)), 
          parnmpc.discrete().t_lift(lift_index), 
          parnmpc.discrete().dtau_lift(lift_index), 
          s[time_stage_before_lift].q, s[time_stage_before_lift].v, 
          s.lift[lift_index], s[time_stage_before_lift+1], 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


double ParNMPCLinearizer::KKTError(const ParNMPC& parnmpc, 
                                   const KKTResidual& kkt_residual) {
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  const int N_all = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_-1) {
      kkt_error_.coeffRef(i) 
          = parnmpc[i].squaredNormKKTResidual(kkt_residual[i], 
                                              parnmpc.discrete().dtau(i));
    }
    else if (i == N_-1) {
      kkt_error_.coeffRef(i) 
          = parnmpc.terminal.squaredNormKKTResidual(
              kkt_residual[i], parnmpc.discrete().dtau(i));
    }
    else if (i < N_+N_impulse) {
      const int impulse_index  = i - N_;
      const int time_stage_before_impulse 
          = parnmpc.discrete().timeStageBeforeImpulse(impulse_index);
      kkt_error_.coeffRef(i) 
          = parnmpc.impulse[impulse_index].squaredNormKKTResidual(
              kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index  = i - (N_+N_impulse);
      kkt_error_.coeffRef(i) 
          = parnmpc.aux[impulse_index].squaredNormKKTResidual(
                kkt_residual.aux[impulse_index], 
                parnmpc.discrete().dtau_aux(impulse_index));
    }
    else {
      const int lift_index = i - (N_+2*N_impulse);
      kkt_error_.coeffRef(i) 
          = parnmpc.lift[lift_index].squaredNormKKTResidual(
              kkt_residual.lift[lift_index], 
              parnmpc.discrete().dtau_lift(lift_index));
    }
  }
  return std::sqrt(kkt_error_.head(N_all).sum());
}


void ParNMPCLinearizer::integrateSolution(ParNMPC& parnmpc, 
                                          const std::vector<Robot>& robots, 
                                          const KKTMatrix& kkt_matrix, 
                                          const KKTResidual& kkt_residual, 
                                          const double primal_step_size, 
                                          const double dual_step_size, 
                                          const Direction& d, 
                                          Solution& s) const {
  assert(robots.size() == nthreads_);
  const int N_impulse = parnmpc.discrete().numImpulseStages();
  const int N_lift = parnmpc.discrete().numLiftStages();
  const int N_all = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_-1) {
      parnmpc[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                              d[i], s[i]);
      parnmpc[i].updateDual(dual_step_size);
    }
    else if (i == N_-1) {
      parnmpc.terminal.updatePrimal(robots[omp_get_thread_num()],  
                                    primal_step_size, d[i], s[i]);
      parnmpc.terminal.updateDual(dual_step_size);
    }
    else if (i < N_+N_impulse) {
      const int impulse_index  = i - N_;
      parnmpc.impulse[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                                  primal_step_size, 
                                                  d.impulse[impulse_index], 
                                                  s.impulse[impulse_index]);
      parnmpc.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index  = i - (N_+N_impulse);
      parnmpc.aux[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                              primal_step_size, 
                                              d.aux[impulse_index], 
                                              s.aux[impulse_index]);
      parnmpc.aux[impulse_index].updateDual(dual_step_size);
    }
    else {
      const int lift_index = i - (N_+2*N_impulse);
      parnmpc.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                            primal_step_size, 
                                            d.lift[lift_index], 
                                            s.lift[lift_index]);
      parnmpc.lift[lift_index].updateDual(dual_step_size);
    }
  }
}

} // namespace idocp