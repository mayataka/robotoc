#include "idocp/ocp/ocp_linearizer.hpp"


namespace idocp {

void OCPLinearizer::linearizeOCP(std::vector<Robot>& robots,
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
    if (i == 0) {
      if (contact_sequence.existImpulseStage(0)) {
        constexpr int impulse_index = 0;
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index);
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps_[0].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(0), 
                                    t, dtau_impulse, q, s[0], s.impulse[impulse_index],
                                    kkt_matrix[0], kkt_residual[0]);
      }
      else if (contact_sequence.existLiftStage(0)) {
        constexpr int lift_index = 0;
        const double dtau_lift = contact_sequence.liftTime(lift_index);
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps_[0].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(0), 
                                    t, dtau_lift, q, s[0], s.lift[lift_index],
                                    kkt_matrix[0], kkt_residual[0]);
      }
      else {
        split_ocps_[0].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(0), 
                                    t, dtau_, q, s[0], s[1],
                                    kkt_matrix[0], kkt_residual[0]);
      }
    }
    else if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps_[i].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(i), 
                                    t+i*dtau_, dtau_impulse, 
                                    q_prev(contact_sequence, s, i),  
                                    s[i], s.impulse[impulse_index],
                                    kkt_matrix[i], kkt_residual[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps_[i].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(i), 
                                    t+i*dtau_, dtau_lift, 
                                    q_prev(contact_sequence, s, i),  
                                    s[i], s.lift[lift_index],
                                    kkt_matrix[i], kkt_residual[i]);
      }
      else {
        split_ocps_[i].linearizeOCP(robots[omp_get_thread_num()], 
                                    contact_sequence.contactStatus(i), 
                                    t+i*dtau_, dtau_, 
                                    q_prev(contact_sequence, s, i),  
                                    s[i], s[i+1],
                                    kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N_) {
      terminal_ocp_.linearizeOCP(robots[omp_get_thread_num()], t+T_, s[N_],
                                 kkt_matrix[N_], kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      split_ocps_.impulse[impulse_index].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.impulseStatus(impulse_index), 
          t+contact_sequence.impulseTime(impulse_index), 
          s[time_stage_before_impulse].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix.impulse[impulse_index], kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      const double dtau_aux 
          = (time_stage_before_impulse+1) * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      split_ocps_.aux[impulse_index].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(time_stage_before_impulse+1), 
          t+contact_sequence.impulseTime(impulse_index), dtau_aux,
          s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[time_stage_before_impulse+1], 
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_before_lift 
          = contact_sequence.timeStageBeforeLift(lift_index);
      const double dtau_aux
          = (time_stage_before_lift+1) * dtau_ 
              - contact_sequence.liftTime(lift_index);
      split_ocps_.lift[lift_index].linearizeOCP(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(time_stage_before_lift+1), 
          t+contact_sequence.liftTime(lift_index), dtau_aux, 
          s[time_stage_before_lift].q, 
          s.lift[lift_index], s[time_stage_before_lift+1],
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


void OCPLinearizer::computeKKTResidual(std::vector<Robot>& robots, 
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
    if (i == 0) {
      if (contact_sequence.existImpulseStage(0)) {
        constexpr int impulse_index = 0;
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index);
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps_[0].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(0), 
                                          t, dtau_impulse, q, s[0], 
                                          s.impulse[impulse_index], 
                                          kkt_matrix[0], kkt_residual[0]);
      }
      else if (contact_sequence.existLiftStage(0)) {
        constexpr int lift_index = 0;
        const double dtau_lift = contact_sequence.liftTime(lift_index);
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps_[0].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(0), 
                                          t, dtau_lift, q, s[0], s.lift[lift_index],
                                          kkt_matrix[0], kkt_residual[0]);
      }
      else {
        split_ocps_[0].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(0), 
                                          t, dtau_, q, s[0], s[1],
                                          kkt_matrix[0], kkt_residual[0]);
      }
    }
    else if (i < N_) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(i), 
                                          t+i*dtau_, dtau_impulse, 
                                          q_prev(contact_sequence, s, i),  
                                          s[i], s.impulse[impulse_index],
                                          kkt_matrix[i], kkt_residual[i]);
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(i), 
                                          t+i*dtau_, dtau_lift, 
                                          q_prev(contact_sequence, s, i),  
                                          s[i], s.lift[lift_index],
                                          kkt_matrix[i], kkt_residual[i]);
      }
      else {
        split_ocps_[i].computeKKTResidual(robots[omp_get_thread_num()], 
                                          contact_sequence.contactStatus(i), 
                                          t+i*dtau_, dtau_, 
                                          q_prev(contact_sequence, s, i),  
                                          s[i], s[i+1],
                                          kkt_matrix[i], kkt_residual[i]);
      }
    }
    else if (i == N_) {
      terminal_ocp_.computeKKTResidual(robots[omp_get_thread_num()], t+T_, s[N_],
                                       kkt_residual[N_]);
    }
    else if (i < N_ + 1 + N_impulse) {
      const int impulse_index  = i - (N_+1);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      split_ocps_.impulse[impulse_index].computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.impulseStatus(impulse_index), 
          t+contact_sequence.impulseTime(impulse_index), 
          s[time_stage_before_impulse].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix.impulse[impulse_index], kkt_residual.impulse[impulse_index]);
    }
    else if (i < N_ + 1 + 2*N_impulse) {
      const int impulse_index  = i - (N_+1+N_impulse);
      const int time_stage_before_impulse 
          = contact_sequence.timeStageBeforeImpulse(impulse_index);
      const double dtau_aux 
          = (time_stage_before_impulse+1) * dtau_ 
              - contact_sequence.impulseTime(impulse_index);
      split_ocps_.aux[impulse_index].computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(time_stage_before_impulse+1), 
          t+contact_sequence.impulseTime(impulse_index), dtau_aux,
          s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[time_stage_before_impulse+1], 
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
    else {
      const int lift_index = i - (N_+1+2*N_impulse);
      const int time_stage_before_lift 
          = contact_sequence.timeStageBeforeLift(lift_index);
      const double dtau_aux
          = (time_stage_before_lift+1) * dtau_ 
              - contact_sequence.liftTime(lift_index);
      split_ocps_.lift[lift_index].computeKKTResidual(
          robots[omp_get_thread_num()], 
          contact_sequence.contactStatus(time_stage_before_lift+1), 
          t+contact_sequence.liftTime(lift_index), dtau_aux, 
          s[time_stage_before_lift].q, 
          s.lift[lift_index], s[time_stage_before_lift+1],
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}

} // namespace idocp 