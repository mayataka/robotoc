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


inline OCPLinearizer::OCPLinearizer()
  : T_(0),
    dtau_(0),
    N_(0),
    num_proc_(0),
    kkt_error_() {
}


inline OCPLinearizer::~OCPLinearizer() {
}


inline void OCPLinearizer::linearizeOCP(
    HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
    std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const HybridSolution& s, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual) const {
  runParallel<internal::LinearizeOCP>(split_ocps, terminal_ocp, robots, 
                                      contact_sequence, t, q, v, s, 
                                      kkt_matrix, kkt_residual);
}


inline void OCPLinearizer::computeKKTResidual(
    HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
    std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
    const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
    const HybridSolution& s, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual) const {
  runParallel<internal::ComputeKKTResidual>(split_ocps, terminal_ocp, robots, 
                                            contact_sequence, t, q, v, s, 
                                            kkt_matrix, kkt_residual);
}


inline double OCPLinearizer::KKTError(const HybridOCP& split_ocps, 
                                      const TerminalOCP& terminal_ocp, 
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
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        kkt_error_.coeffRef(i) = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], dtau_impulse);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        kkt_error_.coeffRef(i) = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], dtau_lift);
      }
      else {
        kkt_error_.coeffRef(i) = split_ocps[i].squaredNormKKTResidual(kkt_residual[i], dtau_);
      }
    }
    else if (i == N_) {
      kkt_error_.coeffRef(N_) = terminal_ocp.squaredNormKKTResidual(kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      kkt_error_.coeffRef(i) 
          = split_ocps.impulse[impulse_index].squaredNormKKTResidual(kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      const double dtau_aux 
          = (time_stage_before_impulse+1) * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      kkt_error_.coeffRef(i) 
          = split_ocps.aux[impulse_index].squaredNormKKTResidual(kkt_residual.aux[impulse_index], dtau_aux);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_before_lift 
          = contact_sequence.timeStageBeforeLift(lift_index);
      const double dtau_aux
          = (time_stage_before_lift+1) * dtau_ 
              - contact_sequence.liftTime(lift_index);
      kkt_error_.coeffRef(i) 
          = split_ocps.lift[lift_index].squaredNormKKTResidual(kkt_residual.lift[lift_index], dtau_aux);
    }
  }
  return std::sqrt(kkt_error_.head(N_all).sum());
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
                                       HybridKKTResidual& kkt_residual) const {
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

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 