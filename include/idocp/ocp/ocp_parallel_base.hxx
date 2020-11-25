#ifndef IDOCP_OCP_PARALLEL_BASE_HXX_ 
#define IDOCP_OCP_PARALLEL_BASE_HXX_

#include "idocp/ocp/ocp_parallel_base.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

namespace idocp {

template <typename... Args>
inline void OCPParallelBase<Derived>::runParallel(Args... args) const {
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
        static_cast<Derived*>(this)->run(i, impulse_index, args);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift 
            = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        static_cast<Derived*>(this)->run(i, lift_index, args);
      }
      else {
        Derived::run(i, args);
      }
    }
    else if (i == N_) {
      Derived::run_terminal(N_, args);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      Derived::run_impulse(impulse_index);
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
      Derived::run_aux(impulse_index);
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
      Derived::run_lift(lift_index);
    }
  }
}


inline const Eigen::VectorXd& OCPParallelBase::q_prev(
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


inline double OCPParallelBase::dtau(const ContactSequence& contact_sequence, 
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


inline bool OCPParallelBase::is_state_constraint_valid(
    const int time_stage_before_impulse) {
  assert(time_stage_before_impulse >= 0);
  if (time_stage_before_impulse > 0) {
    return true
  }
  else {
    return false;
  }
}

} // namespace idocp 

#endif // IDOCP_OCP_LINEARIZER_HXX_ 