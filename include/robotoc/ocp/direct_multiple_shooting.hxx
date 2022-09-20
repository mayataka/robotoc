#ifndef ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HXX_ 
#define ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HXX_

#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include <omp.h>
#include <cassert>


namespace robotoc {
namespace internal {

struct ComputeKKTResidual {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, grid_info, q_prev, 
                                 s, s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, const SplitSolution& s_next, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual,
                         const ImpulseStatus& impulse_status, 
                         const GridInfo& grid_info_next, 
                         SwitchingConstraintJacobian& sc_jacobian,
                         SwitchingConstraintResidual& sc_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, grid_info, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual, impulse_status, 
                                 grid_info_next, sc_jacobian, sc_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev,
                         const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, grid_info, q_prev, s, kkt_matrix, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const ImpulseSplitSolution& s, const SplitSolution& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    impulse_split_ocp.computeKKTResidual(robot, impulse_status, grid_info, 
                                         q_prev, s, s_next, kkt_matrix, kkt_residual);
  }
};


struct ComputeKKTSystem {
  template <typename SplitSolutionType>
  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, const SplitSolutionType& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    split_ocp.computeKKTSystem(robot, contact_status, grid_info, q_prev, 
                               s, s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(SplitOCP& split_ocp, Robot& robot, 
                         const ContactStatus& contact_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, const SplitSolution& s_next, 
                         SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual,
                         const ImpulseStatus& impulse_status,  
                         const GridInfo& grid_info_next, 
                         SwitchingConstraintJacobian& sc_jacobian,
                         SwitchingConstraintResidual& sc_residual) {
    split_ocp.computeKKTSystem(robot, contact_status, grid_info, q_prev, s, 
                               s_next, kkt_matrix, kkt_residual, impulse_status, 
                               grid_info_next, sc_jacobian, sc_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    terminal_ocp.computeKKTSystem(robot, grid_info, q_prev, s, kkt_matrix, kkt_residual);
  }

  static inline void run(ImpulseSplitOCP& impulse_split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, 
                         const GridInfo& grid_info, const Eigen::VectorXd& q_prev, 
                         const ImpulseSplitSolution& s, const SplitSolution& s_next, 
                         SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual) {
    impulse_split_ocp.computeKKTSystem(robot, impulse_status, grid_info, q_prev,  
                                       s, s_next, kkt_matrix, kkt_residual);
  }
};

} // namespace internal
} // namespace robotoc


namespace robotoc {

template <typename Algorithm>
inline void DirectMultipleShooting::runParallel(
    OCP& ocp, aligned_vector<Robot>& robots, 
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_all; ++i) {
    if (i < N) {
      if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i)) {
        assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
        const int contact_phase = ocp.timeDiscretization().contactPhase(i);
        Algorithm::run(ocp[i], robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(contact_phase), 
                       ocp.timeDiscretization().gridInfo(i), q_prev(ocp, q, s, i), s[i], 
                       s.impulse[ocp.timeDiscretization().impulseIndexAfterTimeStage(i)], 
                       kkt_matrix[i], kkt_residual[i]);
        SplitOCP::correctSTOSensitivities(kkt_matrix[i], kkt_residual[i], 
                                          ocp.timeDiscretization().N_phase(contact_phase));
      }
      else if (ocp.timeDiscretization().isTimeStageBeforeLift(i)) {
        assert(!ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1));
        const int contact_phase = ocp.timeDiscretization().contactPhase(i);
        Algorithm::run(ocp[i], robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(contact_phase), 
                       ocp.timeDiscretization().gridInfo(i), q_prev(ocp, q, s, i), s[i], 
                       s.lift[ocp.timeDiscretization().liftIndexAfterTimeStage(i)], 
                       kkt_matrix[i], kkt_residual[i]);
        SplitOCP::correctSTOSensitivities(kkt_matrix[i], kkt_residual[i], 
                                          ocp.timeDiscretization().N_phase(contact_phase));
      }
      else if (ocp.timeDiscretization().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index  
            = ocp.timeDiscretization().impulseIndexAfterTimeStage(i+1);
        const int contact_phase = ocp.timeDiscretization().contactPhase(i);
        Algorithm::run(ocp[i], robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(contact_phase), 
                       ocp.timeDiscretization().gridInfo(i), q_prev(ocp, q, s, i), s[i], 
                       s[i+1], kkt_matrix[i], kkt_residual[i], 
                       contact_sequence->impulseStatus(impulse_index), 
                       ocp.timeDiscretization().gridInfo(i+1), 
                       kkt_matrix.switching[impulse_index], 
                       kkt_residual.switching[impulse_index]);
        SplitOCP::correctSTOSensitivities(kkt_matrix[i], kkt_residual[i], 
                                          kkt_matrix.switching[impulse_index],
                                          ocp.timeDiscretization().N_phase(contact_phase));
      }
      else {
        const int contact_phase = ocp.timeDiscretization().contactPhase(i);
        Algorithm::run(ocp[i], robots[omp_get_thread_num()], 
                       contact_sequence->contactStatus(contact_phase), 
                       ocp.timeDiscretization().gridInfo(i), q_prev(ocp, q, s, i), s[i], 
                       s[i+1], kkt_matrix[i], kkt_residual[i]);
        SplitOCP::correctSTOSensitivities(kkt_matrix[i], kkt_residual[i], 
                                          ocp.timeDiscretization().N_phase(contact_phase));
      }
    }
    else if (i == N) {
      Algorithm::run(ocp.terminal, robots[omp_get_thread_num()], 
                     ocp.timeDiscretization().gridInfo(N), q_prev(ocp, q, s, N), s[N], 
                     kkt_matrix[N], kkt_residual[N]);
    }
    else if (i < N+1+N_impulse) {
      const int impulse_index  = i - (N+1);
      const int time_stage_before_impulse 
          = ocp.timeDiscretization().timeStageBeforeImpulse(impulse_index);
      Algorithm::run(ocp.impulse[impulse_index], robots[omp_get_thread_num()], 
                     contact_sequence->impulseStatus(impulse_index), 
                     ocp.timeDiscretization().gridInfoImpulse(impulse_index), 
                     s[time_stage_before_impulse].q, s.impulse[impulse_index], 
                     s.aux[impulse_index], kkt_matrix.impulse[impulse_index], 
                     kkt_residual.impulse[impulse_index]);
    }
    else if (i < N+1+2*N_impulse) {
      const int impulse_index  = i - (N+1+N_impulse);
      const int time_stage_after_impulse 
          = ocp.timeDiscretization().timeStageAfterImpulse(impulse_index);
      const int contact_phase = ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index);
      Algorithm::run(ocp.aux[impulse_index], robots[omp_get_thread_num()], 
                     contact_sequence->contactStatus(contact_phase),
                     ocp.timeDiscretization().gridInfoAux(impulse_index), s.impulse[impulse_index].q, 
                     s.aux[impulse_index], s[time_stage_after_impulse], 
                     kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
      SplitOCP::correctSTOSensitivities(kkt_matrix.aux[impulse_index], 
                                        kkt_residual.aux[impulse_index], 
                                        ocp.timeDiscretization().N_phase(contact_phase));
    }
    else {
      const int lift_index = i - (N+1+2*N_impulse);
      const int time_stage_after_lift
          = ocp.timeDiscretization().timeStageAfterLift(lift_index);
      const int contact_phase = ocp.timeDiscretization().contactPhaseAfterLift(lift_index);
      Algorithm::run(ocp.lift[lift_index], robots[omp_get_thread_num()], 
                     contact_sequence->contactStatus(contact_phase),
                     ocp.timeDiscretization().gridInfoLift(lift_index), s[time_stage_after_lift-1].q, 
                     s.lift[lift_index], s[time_stage_after_lift], 
                     kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
      SplitOCP::correctSTOSensitivities(kkt_matrix.lift[lift_index], 
                                        kkt_residual.lift[lift_index],
                                        ocp.timeDiscretization().N_phase(contact_phase));
    }
  }
}


inline const Eigen::VectorXd& DirectMultipleShooting::q_prev(const OCP& ocp, 
                                                             const Eigen::VectorXd& q, 
                                                             const Solution& s, 
                                                             const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage <= ocp.timeDiscretization().N());
  if (time_stage == 0) {
    return q;
  }
  else if (ocp.timeDiscretization().isTimeStageBeforeImpulse(time_stage-1)) {
    return s.aux[ocp.timeDiscretization().impulseIndexAfterTimeStage(time_stage-1)].q;
  }
  else if (ocp.timeDiscretization().isTimeStageBeforeLift(time_stage-1)) {
    return s.lift[ocp.timeDiscretization().liftIndexAfterTimeStage(time_stage-1)].q;
  }
  else {
    return s[time_stage-1].q;
  }
}


inline double DirectMultipleShooting::dts_stage(const OCP& ocp, 
                                                const Direction& d, 
                                                const int time_stage) {
  const int phase = ocp.timeDiscretization().contactPhase(time_stage);
  return ((d[time_stage].dts_next-d[time_stage].dts) 
            / ocp.timeDiscretization().N_phase(phase));
}


inline double DirectMultipleShooting::dts_aux(const OCP& ocp, 
                                              const Direction& d, 
                                              const int impulse_index) {
  const int phase = ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index);
  return ((d.aux[impulse_index].dts_next-d.aux[impulse_index].dts) 
            / ocp.timeDiscretization().N_phase(phase));
}


inline double DirectMultipleShooting::dts_lift(const OCP& ocp, 
                                               const Direction& d, 
                                               const int lift_index) {
  const int phase = ocp.timeDiscretization().contactPhaseAfterLift(lift_index);
  return ((d.lift[lift_index].dts_next-d.lift[lift_index].dts) 
            / ocp.timeDiscretization().N_phase(phase));
}

} // namespace robotoc 

#endif // ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HXX_ 