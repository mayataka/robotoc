#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp{

OCPLinearizer::OCPLinearizer(const double T, const int N, 
                             const int max_num_impulse, 
                             const int num_proc) 
  : T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    kkt_error_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


OCPLinearizer::OCPLinearizer()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0),
    kkt_error_() {
}


OCPLinearizer::~OCPLinearizer() {
}


void OCPLinearizer::initConstraints(HybridOCP& split_ocps, 
                                    std::vector<Robot>& robots, 
                                    const ContactSequence& contact_sequence, 
                                    const HybridSolution& s) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse 
            = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps[i].initConstraints(robots[omp_get_thread_num()], i, 
                                      dtau_impulse, s[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift 
            = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps[i].initConstraints(robots[omp_get_thread_num()], i, 
                                      dtau_lift, s[i]);
      }
      else {
        split_ocps[i].initConstraints(robots[omp_get_thread_num()], i, 
                                      dtau_, s[i]);
      }
    }
    else if (i == N_) {
      split_ocps.terminal.initConstraints(robots[omp_get_thread_num()], N_, 
                                          dtau_, s[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      split_ocps.impulse[impulse_index].initConstraints(
          robots[omp_get_thread_num()], s.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = contact_sequence.timeStageAfterImpulse(impulse_index);
      const double dtau_aux 
          = time_stage_after_impulse * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      split_ocps.aux[impulse_index].initConstraints(
          robots[omp_get_thread_num()], 0, dtau_aux, s.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift 
          = contact_sequence.timeStageAfterLift(lift_index);
      const double dtau_aux
          = time_stage_after_lift * dtau_ 
              - contact_sequence.liftTime(lift_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      split_ocps.lift[lift_index].initConstraints(
          robots[omp_get_thread_num()], 0, dtau_aux, s.lift[lift_index]);
    }
  }
}


void OCPLinearizer::linearizeOCP(HybridOCP& split_ocps, 
                                 std::vector<Robot>& robots, 
                                 const ContactSequence& contact_sequence, 
                                 const double t, const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const HybridSolution& s, 
                                 HybridKKTMatrix& kkt_matrix, 
                                 HybridKKTResidual& kkt_residual) const {
  runParallel<internal::LinearizeOCP>(split_ocps, robots, contact_sequence, 
                                      t, q, v, s, kkt_matrix, kkt_residual);
}


void OCPLinearizer::computeKKTResidual(HybridOCP& split_ocps, 
                                       std::vector<Robot>& robots, 
                                       const ContactSequence& contact_sequence, 
                                       const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const HybridSolution& s, 
                                       HybridKKTMatrix& kkt_matrix, 
                                       HybridKKTResidual& kkt_residual) const {
  runParallel<internal::ComputeKKTResidual>(split_ocps, robots, contact_sequence,
                                            t, q, v, s, kkt_matrix, kkt_residual);
}


double OCPLinearizer::KKTError(const HybridOCP& split_ocps, 
                               const ContactSequence& contact_sequence, 
                               const HybridKKTResidual& kkt_residual) {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse 
            = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        kkt_error_.coeffRef(i) 
            = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], 
                                                   dtau_impulse);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift 
            = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        kkt_error_.coeffRef(i) 
            = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], dtau_lift);
      }
      else {
        kkt_error_.coeffRef(i) 
            = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], dtau_);
      }
    }
    else if (i == N_) {
      kkt_error_.coeffRef(N_) 
          = split_ocps.terminal.squaredNormKKTResidual(kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      kkt_error_.coeffRef(i) 
          = split_ocps.impulse[impulse_index].squaredNormKKTResidual(
              kkt_residual.impulse[impulse_index], 
              is_state_constraint_valid(time_stage_before_impulse));
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = contact_sequence.timeStageAfterImpulse(impulse_index);
      const double dtau_aux 
          = time_stage_after_impulse * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      kkt_error_.coeffRef(i) 
          = split_ocps.aux[impulse_index].squaredNormKKTResidual(
                kkt_residual.aux[impulse_index], dtau_aux);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift 
          = contact_sequence.timeStageAfterLift(lift_index);
      const double dtau_aux
          = time_stage_after_lift * dtau_ 
              - contact_sequence.liftTime(lift_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      kkt_error_.coeffRef(i) 
          = split_ocps.lift[lift_index].squaredNormKKTResidual(
              kkt_residual.lift[lift_index], dtau_aux);
    }
  }
  return std::sqrt(kkt_error_.head(N_all).sum());
}


void OCPLinearizer::integrateSolution(HybridOCP& split_ocps, 
                                      const std::vector<Robot>& robots, 
                                      const ContactSequence& contact_sequence, 
                                      const HybridKKTMatrix& kkt_matrix, 
                                      const HybridKKTResidual& kkt_residual, 
                                      const double primal_step_size, 
                                      const double dual_step_size, 
                                      HybridDirection& d, 
                                      HybridSolution& s) const {
  assert(robots.size() == num_proc_);
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps[i].computeCondensedDualDirection(robots[omp_get_thread_num()], 
                                                    dtau_impulse, kkt_matrix[i], 
                                                    kkt_residual[i], 
                                                    d.impulse[impulse_index], d[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps[i].computeCondensedDualDirection(robots[omp_get_thread_num()], 
                                                    dtau_lift, kkt_matrix[i], 
                                                    kkt_residual[i], 
                                                    d.lift[lift_index], d[i]);
      }
      else {
        split_ocps[i].computeCondensedDualDirection(robots[omp_get_thread_num()], 
                                                    dtau_, kkt_matrix[i], 
                                                    kkt_residual[i], d[i+1], d[i]);
      }
      split_ocps[i].updatePrimal(robots[omp_get_thread_num()], primal_step_size, 
                                 dtau_, d[i], s[i]);
      split_ocps[i].updateDual(dual_step_size);
    }
    else if (i == N_) {
      split_ocps.terminal.updatePrimal(robots[omp_get_thread_num()], 
                                       primal_step_size, d[N_], s[N_]);
      split_ocps.terminal.updateDual(dual_step_size);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      split_ocps.impulse[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index], d.aux[impulse_index], 
          d.impulse[impulse_index]);
      split_ocps.impulse[impulse_index].updatePrimal(
          robots[omp_get_thread_num()], primal_step_size, 
          d.impulse[impulse_index], s.impulse[impulse_index], 
          is_state_constraint_valid(contact_sequence.timeStageBeforeImpulse(impulse_index)));
      split_ocps.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = contact_sequence.timeStageAfterImpulse(impulse_index);
      const double dtau_aux 
          = time_stage_after_impulse * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      split_ocps.aux[impulse_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], dtau_aux,
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          d[time_stage_after_impulse], d.aux[impulse_index]);
      split_ocps.aux[impulse_index].updatePrimal(robots[omp_get_thread_num()], 
                                                 primal_step_size, dtau_aux,
                                                 d.aux[impulse_index], 
                                                 s.aux[impulse_index]);
      split_ocps.aux[impulse_index].updateDual(dual_step_size);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift 
          = contact_sequence.timeStageAfterLift(lift_index);
      const double dtau_aux
          = time_stage_after_lift * dtau_ 
              - contact_sequence.liftTime(lift_index);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      split_ocps.lift[lift_index].computeCondensedDualDirection(
          robots[omp_get_thread_num()], dtau_aux,
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index],
          d[time_stage_after_lift], d.lift[lift_index]);
      split_ocps.lift[lift_index].updatePrimal(robots[omp_get_thread_num()], 
                                               primal_step_size, dtau_aux,
                                               d.lift[lift_index], 
                                               s.lift[lift_index]);
      split_ocps.lift[lift_index].updateDual(dual_step_size);
    }
  }
}

} // namespace idocp