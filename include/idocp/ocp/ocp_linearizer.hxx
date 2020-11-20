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
                         KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
    split_ocp.linearizeOCP(robot, contact_status, t, dtau, q_prev, s, s_next,
                           kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot, 
                         const double t, const SplitSolution& s, 
                         KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
    terminal_ocp.linearizeOCP(robot, t, s, kkt_matrix, kkt_residual);
  }

  template <typename ConfigVectorType>
  static inline void run(SplitImpulseOCP& split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseKKTMatrix& kkt_matrix, 
                         ImpulseKKTResidual& kkt_residual) {
    split_ocp.linearizeOCP(robot, impulse_status, t, q_prev, s, s_next,
                           kkt_matrix, kkt_residual);
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
                         KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, 
                                 s_next, kkt_matrix, kkt_residual);
  }

  static inline void run(TerminalOCP& terminal_ocp, Robot& robot,  
                         const double t, const SplitSolution& s, 
                         KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
    terminal_ocp.computeKKTResidual(robot, t, s, kkt_residual);
  }

  template <typename ConfigVectorType>
  static inline void run(SplitImpulseOCP& split_ocp, Robot& robot, 
                         const ImpulseStatus& impulse_status, const double t, 
                         const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
                         const ImpulseSplitSolution& s, 
                         const SplitSolution& s_next, 
                         ImpulseKKTMatrix& kkt_matrix, 
                         ImpulseKKTResidual& kkt_residual) {
    split_ocp.computeKKTResidual(robot, impulse_status, t, q_prev, s, s_next,
                                 kkt_matrix, kkt_residual);
  }
};

} // namespace internal
} // namespace idocp


namespace idocp {

inline OCPLinearizer::OCPLinearizer(const double T, const int N, 
                                    const int max_num_impulse, 
                                    const int num_proc) 
  : T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc) {
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


inline OCPLinearizer::OCPLinearizer()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0) {
}


inline OCPLinearizer::~OCPLinearizer() {
}


inline void OCPLinearizer::linearizeOCP(
    HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
    std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const HybridSolution& s, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual) {
  runParallel<internal::LinearizeOCP>(split_ocps, terminal_ocp, robots, 
                                      contact_sequence, t, q, v, s, 
                                      kkt_matrix, kkt_residual);
}


inline void OCPLinearizer::computeKKTResidual(
    HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
    std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const HybridSolution& s, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual) {
  runParallel<internal::ComputeKKTResidual>(split_ocps, terminal_ocp, robots, 
                                            contact_sequence, t, q, v, s, 
                                            kkt_matrix, kkt_residual);
}


template <typename Algorithm>
inline void OCPLinearizer::runParallel(HybridOCP& split_ocps, 
                                       TerminalOCP& terminal_ocp,
                                       std::vector<Robot>& robots,
                                       const ContactSequence& contact_sequence,
                                       const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const HybridSolution& s,
                                       HybridKKTMatrix& kkt_matrix,
                                       HybridKKTResidual& kkt_residual) {
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
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        Algorithm::run(split_ocps[i], robots[omp_get_thread_num()], 
                       contact_sequence.contactStatus(i), t+i*dtau_, 
                       dtau_impulse, q_prev(contact_sequence, q, s, i), s[i], 
                       s.impulse[impulse_index], kkt_matrix[i], kkt_residual[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
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
      Algorithm::run(terminal_ocp, robots[omp_get_thread_num()], t+T_, s[N_], 
                     kkt_matrix[N_], kkt_residual[N_]);
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
                     kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      const double dtau_aux 
          = (time_stage_before_impulse+1) * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      Algorithm::run(split_ocps.aux[impulse_index], 
                     robots[omp_get_thread_num()], 
                     contact_sequence.contactStatus(time_stage_before_impulse+1), 
                     t+contact_sequence.impulseTime(impulse_index), dtau_aux,
                     s.impulse[impulse_index].q, s.aux[impulse_index], 
                     s[time_stage_before_impulse+1], 
                     kkt_matrix.aux[impulse_index], 
                     kkt_residual.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_before_lift 
          = contact_sequence.timeStageBeforeLift(lift_index);
      const double dtau_aux
          = (time_stage_before_lift+1) * dtau_ 
              - contact_sequence.liftTime(lift_index);
      Algorithm::run(split_ocps.lift[lift_index], robots[omp_get_thread_num()], 
                     contact_sequence.contactStatus(time_stage_before_lift+1), 
                     t+contact_sequence.liftTime(lift_index), dtau_aux, 
                     s[time_stage_before_lift].q, s.lift[lift_index], 
                     s[time_stage_before_lift+1], 
                     kkt_matrix.lift[lift_index], 
                     kkt_residual.lift[lift_index]);
    }
  }
}


inline const Eigen::VectorXd& OCPLinearizer::q_prev(
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q,
    const HybridSolution& s, const int time_stage) const {
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


// inline void OCPLinearizer::computePrimalDirectionInitial(
//     const RiccatiFactorizer factorizer, 
//     const RiccatiFactorization factorization, SplitDirection& d, 
//     const bool exist_state_constraint) {
//   RiccatiFactorizer::computeCostateDirection(factorization, d, 
//                                              exist_state_constraint);
//   factorizer.computeControlInputDirection(factorization, d, 
//                                           exist_state_constraint);
// }


// template <typename VectorType>
// inline void OCPLinearizer::computePrimalDirection(
//     const RiccatiFactorizer factorizer, const RiccatiFactorization factorization, 
//     const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d, 
//     const bool exist_state_constraint) {
//   RiccatiFactorizer::computeStateDirection(factorization, dx0, d,
//                                            exist_state_constraint);
//   RiccatiFactorizer::computeCostateDirection(factorization, d, 
//                                              exist_state_constraint);
//   factorizer.computeControlInputDirection(factorization, d, 
//                                           exist_state_constraint);
// }


// template <typename VectorType>
// inline void OCPLinearizer::computePrimalDirectionTerminal(
//     const RiccatiFactorization factorization, 
//     const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d) {
//   RiccatiFactorizer::computeStateDirection(factorization, dx0, d, false);
//   RiccatiFactorizer::computeCostateDirection(factorization, d, false);
// }
 

// template <typename VectorType>
// inline void OCPLinearizer::computePrimalDirectionImpulse(
//     const RiccatiFactorization factorization, 
//     const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d) {
//   ImpulseRiccatiFactorizer::computeStateDirection(factorization, dx0, d);
//   ImpulseRiccatiFactorizer::computeCostateDirection(factorization, d);
// }
 

// inline void OCPLinearizer::aggregateLagrangeMultiplierDirection(
//     const ContactSequence& contact_sequence,
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
//     const std::vector<ImpulseSplitDirection>& d_impulse, const int time_stage,
//     RiccatiFactorization& riccati_factorization) {
//   assert(time_stage >= 0);
//   const int num_impulse = contact_sequence.totalNumImpulseStages();
//   riccati_factorization.n.setZero();
//   for (int i=num_impulse-1; i>=0; --i) {
//     if (contact_sequence.timeStageBeforeImpulse(i) < time_stage) {
//       break;
//     }
//     else {
//       riccati_factorization.n.noalias() 
//           += constraint_factorization[i].T(time_stage) * d_impulse[i].dxi();
//     }
//   }
// }


// inline void OCPLinearizer::aggregateLagrangeMultiplierDirectionImpulse(
//     const ContactSequence& contact_sequence,
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
//     const std::vector<ImpulseSplitDirection>& d_impulse, 
//     const int impulse_index,
//     RiccatiFactorization& impulse_riccati_factorization) {
//   assert(impulse_index >= 0);
//   assert(impulse_index < contact_sequence.totalNumImpulseStages());
//   const int num_impulse = contact_sequence.totalNumImpulseStages();
//   impulse_riccati_factorization.n.setZero();
//   for (int i=num_impulse-1; i>=impulse_index; --i) {
//     impulse_riccati_factorization.n.noalias() 
//         += constraint_factorization[i].T_impulse(impulse_index) * d_impulse[i].dxi();
//   }
// }


// inline void OCPLinearizer::aggregateLagrangeMultiplierDirectionAux(
//     const ContactSequence& contact_sequence,
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
//     const std::vector<ImpulseSplitDirection>& d_impulse, 
//     const int impulse_index,
//     RiccatiFactorization& aux_riccati_factorization) {
//   assert(impulse_index >= 0);
//   assert(impulse_index < contact_sequence.totalNumImpulseStages());
//   const int num_impulse = contact_sequence.totalNumImpulseStages();
//   aux_riccati_factorization.n.setZero();
//   for (int i=num_impulse-1; i>=impulse_index; --i) {
//     aux_riccati_factorization.n.noalias() 
//         += constraint_factorization[i].T_aux(impulse_index) * d_impulse[i].dxi();
//   }
// }


// inline void OCPLinearizer::aggregateLagrangeMultiplierDirectionLift(
//     const ContactSequence& contact_sequence,
//     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
//     const std::vector<ImpulseSplitDirection>& d_impulse, 
//     const int lift_index,
//     RiccatiFactorization& lift_riccati_factorization) {
//   assert(lift_index >= 0);
//   assert(lift_index < contact_sequence.totalNumLiftStages());
//   const int num_impulse = contact_sequence.totalNumImpulseStages();
//   const int time_stage_before_lift = contact_sequence.timeStageBeforeLift(lift_index);
//   lift_riccati_factorization.n.setZero();
//   for (int i=num_impulse-1; i>=0; --i) {
//     if (contact_sequence.timeStageBeforeImpulse(i) < time_stage_before_lift) {
//       break;
//     }
//     else {
//       lift_riccati_factorization.n.noalias() 
//           += constraint_factorization[i].T_lift(lift_index) * d_impulse[i].dxi();
//     }
//   }
// }

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 