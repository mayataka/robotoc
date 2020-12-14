#ifndef IDOCP_OCP_LINEARIZER_HXX_ 
#define IDOCP_OCP_LINEARIZER_HXX_

#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {
namespace internal {

struct LinearizeOCP {
  template <typename ConfigVectorType, typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dtau, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const SplitSolution& s, const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    split_ocp.linearizeOCP(robot, contact_status, t, dtau, q_prev, s, s_next,
                           kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot, 
                         const double t, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    terminal_ocp.linearizeOCP(robot, t, s, kkt_matrix, kkt_residual);
  }

  template <typename ConfigVectorType>
  static inline void run(SplitImpulseOCP& split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual,
                         const bool _is_state_constraint_valid) {
    split_ocp.linearizeOCP(robot, impulse_status, t, q_prev, s, s_next,
                           kkt_matrix, kkt_residual, 
                           _is_state_constraint_valid);
  }
};


struct ComputeKKTResidual {
  template <typename ConfigVectorType, typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dtau, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const double t, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, t, s, kkt_residual);
  }

  template <typename ConfigVectorType>
  static inline void run(SplitImpulseOCP& split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual,
                         const bool _is_state_constraint_valid) {
    split_ocp.computeKKTResidual(robot, impulse_status, t, q_prev, s, s_next,
                                 kkt_matrix, kkt_residual, 
                                 _is_state_constraint_valid);
  }
};


struct initConstraints {
  template <typename ConfigVectorType, typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dtau, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const double t, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, t, s, kkt_residual);
  }

  template <typename ConfigVectorType>
  static inline void run(SplitImpulseOCP& split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual,
                         const bool is_state_constraint_valid) {
    split_ocp.computeKKTResidual(robot, impulse_status, t, q_prev, s, s_next,
                                 kkt_matrix, kkt_residual, 
                                 is_state_constraint_valid);
  }
};

} // namespace internal
} // namespace idocp


namespace idocp {

template <typename Algorithm>
inline void OCPLinearizer::runParallel(HybridOCP& split_ocps, 
                                       std::vector<Robot>& robots,
                                       const ContactSequence& contact_sequence,
                                       const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, KKTMatrix& kkt_matrix,
                                       KKTResidual& kkt_residual) const {
  assert(robots.size() == num_proc_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
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
        Algorithm::run(split_ocps[i], robots[omp_get_thread_num()], 
                       contact_sequence.contactStatus(i), t+i*dtau_, 
                       dtau_impulse, q_prev(contact_sequence, q, s, i), s[i], 
                       s.impulse[impulse_index], kkt_matrix[i], kkt_residual[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift 
            = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        Algorithm::run(split_ocps[i], robots[omp_get_thread_num()], 
                       contact_sequence.contactStatus(i), t+i*dtau_, 
                       dtau_lift, q_prev(contact_sequence, q, s, i), s[i], 
                       s.lift[lift_index], kkt_matrix[i], kkt_residual[i]);
      }
      else {
        Algorithm::run(split_ocps[i], robots[omp_get_thread_num()], 
                       contact_sequence.contactStatus(i), t+i*dtau_, 
                       dtau_, q_prev(contact_sequence, q, s, i), s[i], s[i+1], 
                       kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N_) {
      Algorithm::run(split_ocps.terminal, robots[omp_get_thread_num()], t+T_, 
                     s[N_], kkt_matrix[N_], kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      Algorithm::run(split_ocps.impulse[impulse_index],
                     robots[omp_get_thread_num()], 
                     contact_sequence.impulseStatus(impulse_index), 
                     t+contact_sequence.impulseTime(impulse_index), 
                     s[time_stage_before_impulse].q, s.impulse[impulse_index], 
                     s.aux[impulse_index], 
                     kkt_matrix.impulse[impulse_index], 
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
      Algorithm::run(split_ocps.aux[impulse_index], 
                     robots[omp_get_thread_num()], 
                     contact_sequence.contactStatus(time_stage_after_impulse), 
                     t+contact_sequence.impulseTime(impulse_index), dtau_aux,
                     s.impulse[impulse_index].q, s.aux[impulse_index], 
                     s[time_stage_after_impulse], 
                     kkt_matrix.aux[impulse_index], 
                     kkt_residual.aux[impulse_index]);
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
      Algorithm::run(split_ocps.lift[lift_index], robots[omp_get_thread_num()], 
                     contact_sequence.contactStatus(time_stage_after_lift), 
                     t+contact_sequence.liftTime(lift_index), dtau_aux, 
                     s[time_stage_after_lift-1].q, s.lift[lift_index], 
                     s[time_stage_after_lift], 
                     kkt_matrix.lift[lift_index], 
                     kkt_residual.lift[lift_index]);
    }
  }
}


inline const Eigen::VectorXd& OCPLinearizer::q_prev(
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q,
    const Solution& s, const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage <= N_);
  if (time_stage == 0) {
    return q;
  }
  else if (contact_sequence.existImpulseStage(time_stage-1)) {
    return s.aux[contact_sequence.impulseIndex(time_stage-1)].q;
  }
  else if (contact_sequence.existLiftStage(time_stage-1)) {
    return s.lift[contact_sequence.liftIndex(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}


inline double OCPLinearizer::dtau(const ContactSequence& contact_sequence, 
                                  const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  if (contact_sequence.existImpulseStage(time_stage)) {
    return (contact_sequence.impulseTime(contact_sequence.impulseIndex(time_stage)) 
              - time_stage * dtau_);
  }
  else if (contact_sequence.existLiftStage(time_stage)) {
    return (contact_sequence.liftTime(contact_sequence.liftIndex(time_stage)) 
              - time_stage * dtau_);
  }
  else {
    return dtau_;
  }
}


inline bool OCPLinearizer::is_state_constraint_valid(
    const int time_stage_before_impulse) {
  assert(time_stage_before_impulse >= 0);
  if (time_stage_before_impulse > 0) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 