#ifndef IDOCP_OCP_LINEARIZER_HXX_ 
#define IDOCP_OCP_LINEARIZER_HXX_

#include "idocp/ocp/ocp_linearizer.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {
namespace internal {

struct LinearizeOCP {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dtau, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.linearizeOCP(robot, contact_status, t, dtau, q_prev, s, s_next,
                           kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot, 
                         const double t, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.linearizeOCP(robot, t, s, kkt_matrix, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::VectorXd& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual,
                         const bool is_state_constraint_valid) {
    impulse_split_ocp.linearizeOCP(robot, impulse_status, t, q_prev, s, s_next,
                                   kkt_matrix, kkt_residual, 
                                   is_state_constraint_valid);
  }
};


struct ComputeKKTResidual {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, const double t, 
                         const double dtau, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, 
                         const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const double t, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, t, s, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::VectorXd& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual,
                         const bool is_state_constraint_valid) {
    impulse_split_ocp.computeKKTResidual(robot, impulse_status, t, q_prev, s, 
                                         s_next, kkt_matrix, kkt_residual, 
                                         is_state_constraint_valid);
  }
};

} // namespace internal
} // namespace idocp


namespace idocp {

template <typename Algorithm>
inline void OCPLinearizer::runParallel(OCP& ocp, 
                                       const OCPDiscretizer& ocp_discretizer, 
                                       std::vector<Robot>& robots,
                                       const ContactSequence& contact_sequence,
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, KKTMatrix& kkt_matrix,
                                       KKTResidual& kkt_residual) const {
  assert(robots.size() == num_proc_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N_impulse = ocp_discretizer.numImpulseStages();
  const int N_lift = ocp_discretizer.numLiftStages();
  const int N_all = N_ + 1 + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)), 
            ocp_discretizer.t(i), ocp_discretizer.dtau(i), 
            q_prev(ocp_discretizer, q, s, i), 
            s[i], s.impulse[ocp_discretizer.impulseIndex(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)), 
            ocp_discretizer.t(i), ocp_discretizer.dtau(i), 
            q_prev(ocp_discretizer, q, s, i), 
            s[i], s.lift[ocp_discretizer.liftIndex(i)], 
            kkt_matrix[i], kkt_residual[i]);
      }
      else {
        Algorithm::run(
            ocp[i], robots[omp_get_thread_num()], 
            contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)), 
            ocp_discretizer.t(i), ocp_discretizer.dtau(i), 
            q_prev(ocp_discretizer, q, s, i), 
            s[i], s[i+1], kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N_) {
      Algorithm::run(ocp.terminal, robots[omp_get_thread_num()], 
                     ocp_discretizer.t(N_), s[N_], 
                     kkt_matrix[N_], kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = ocp_discretizer.timeStageBeforeImpulse(impulse_index);
      const bool is_state_constraint_valid = (time_stage_before_impulse > 0);
      Algorithm::run(ocp.impulse[impulse_index], robots[omp_get_thread_num()], 
                     contact_sequence.impulseStatus(impulse_index), 
                     ocp_discretizer.t_impulse(impulse_index), 
                     s[time_stage_before_impulse].q, s.impulse[impulse_index], 
                     s.aux[impulse_index], kkt_matrix.impulse[impulse_index], 
                     kkt_residual.impulse[impulse_index], 
                     is_state_constraint_valid);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp_discretizer.timeStageAfterImpulse(impulse_index);
      Algorithm::run(
          ocp.aux[impulse_index], robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              ocp_discretizer.contactPhaseAfterImpulse(impulse_index)), 
          ocp_discretizer.t_impulse(impulse_index), 
          ocp_discretizer.dtau_aux(impulse_index), s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[time_stage_after_impulse], 
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp_discretizer.timeStageAfterLift(lift_index);
      Algorithm::run(
          ocp.lift[lift_index], robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(
              ocp_discretizer.contactPhaseAfterLift(lift_index)), 
          ocp_discretizer.t_lift(lift_index), 
          ocp_discretizer.dtau_lift(lift_index), s[time_stage_after_lift-1].q, 
          s.lift[lift_index], s[time_stage_after_lift], 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


inline const Eigen::VectorXd& OCPLinearizer::q_prev(
    const OCPDiscretizer& ocp_discretizer, const Eigen::VectorXd& q, 
    const Solution& s, const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage <= ocp_discretizer.N());
  if (time_stage == 0) {
    return q;
  }
  else if (ocp_discretizer.isTimeStageBeforeImpulse(time_stage-1)) {
    return s.aux[ocp_discretizer.impulseIndex(time_stage-1)].q;
  }
  else if (ocp_discretizer.isTimeStageBeforeLift(time_stage-1)) {
    return s.lift[ocp_discretizer.liftIndex(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 