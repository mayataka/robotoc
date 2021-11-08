#ifndef ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HXX_
#define ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HXX_

#include "robotoc/hybrid/switching_time_optimization.hpp"

#include <cassert>


namespace robotoc {

inline SwitchingTimeOptimization::SwitchingTimeOptimization(
    const std::shared_ptr<STOCostFunction>& sto_cost, 
    const std::shared_ptr<STOConstraints>& sto_constraints, 
    const int max_num_impulse_events)
  : sto_cost_(sto_cost), 
    sto_constraints_(sto_constraints),
    sto_reg_(STORegularization::defaultSTORegularization()),
    max_num_impulse_events_(max_num_impulse_events),
    kkt_error_(0),
    cost_val_(0),
    h_phase_(Eigen::VectorXd::Zero(2*max_num_impulse_events+1)),
    is_sto_enabled_(true) {
}


inline SwitchingTimeOptimization::SwitchingTimeOptimization() 
  : sto_cost_(), 
    sto_constraints_(),
    sto_reg_(),
    max_num_impulse_events_(0),
    kkt_error_(0),
    cost_val_(0),
    h_phase_(),
    is_sto_enabled_(false) {
}


inline void SwitchingTimeOptimization::initConstraints(const OCP& ocp) const {
  if (is_sto_enabled_) {
    sto_constraints_->setSlack(ocp.discrete());
  }
}


inline void SwitchingTimeOptimization::computeKKTResidual(
    const OCP& ocp, KKTResidual& kkt_residual) {
  if (is_sto_enabled_) {
    cost_val_ = sto_cost_->linearizeCost(ocp.discrete(), kkt_residual); 
    sto_constraints_->linearizeConstraints(ocp.discrete(), kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
  }
}


inline void SwitchingTimeOptimization::computeKKTSystem(
    const OCP& ocp, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  if (is_sto_enabled_) {
    cost_val_ = sto_cost_->quadratizeCost(ocp.discrete(), kkt_matrix, 
                                          kkt_residual); 
    sto_constraints_->linearizeConstraints(ocp.discrete(), kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
    sto_constraints_->condenseSlackAndDual(ocp.discrete(), kkt_matrix, 
                                           kkt_residual);
  }
}


inline void SwitchingTimeOptimization::applyRegularization(
    const OCP& ocp, const double kkt_error, KKTMatrix& kkt_matrix) const {
  if (is_sto_enabled_) {
    sto_reg_.applyRegularization(ocp, kkt_error, kkt_matrix);
  }
}


inline void SwitchingTimeOptimization::computeDirection(const OCP& ocp, 
                                                        const Direction& d) {
  if (is_sto_enabled_) {
    sto_constraints_->expandSlackAndDual(ocp.discrete(), d);
  }
}


inline double SwitchingTimeOptimization::maxPrimalStepSize() const {
  if (is_sto_enabled_) {
    return sto_constraints_->maxPrimalStepSize();
  }
  else {
    return 1.0;
  }
}


inline double SwitchingTimeOptimization::maxDualStepSize() const {
  if (is_sto_enabled_) {
    return sto_constraints_->maxDualStepSize();
  }
  else {
    return 1.0;
  }
}


inline double SwitchingTimeOptimization::KKTError() const {
  return kkt_error_;
}


inline double SwitchingTimeOptimization::KKTError(
    const OCP& ocp, const KKTResidual& kkt_residual) {
  const int N = ocp.discrete().N();
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  const int N_all = N + 1 + 3*N_impulse + N_lift;
  h_phase_.setZero();
  for (int stage=0; stage<N; ++stage) {
    h_phase_.coeffRef(ocp.discrete().contactPhase(stage))
        += kkt_residual[stage].h;
  }
  for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
    h_phase_.coeffRef(ocp.discrete().contactPhaseAfterImpulse(impulse_index))
        += kkt_residual.aux[impulse_index].h;
  }
  for (int lift_index=0; lift_index<N_lift; ++lift_index) {
    h_phase_.coeffRef(ocp.discrete().contactPhaseAfterLift(lift_index))
        += kkt_residual.lift[lift_index].h;
  }
  double kkt_error = 0;
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<N_impulse+N_lift; ++event_index) {
    if (ocp.discrete().eventType(event_index) == DiscreteEventType::Impulse) {
      if (ocp.discrete().isSTOEnabledImpulse(impulse_index)) {
        const double hdiff = h_phase_.coeff(event_index) 
                              - h_phase_.coeff(event_index+1);
        kkt_error += hdiff * hdiff;
      }
      ++impulse_index;
    }
    else {
      assert(ocp.discrete().eventType(event_index) == DiscreteEventType::Lift);
      if (ocp.discrete().isSTOEnabledLift(lift_index)) {
        const double hdiff = h_phase_.coeff(event_index) 
                              - h_phase_.coeff(event_index+1);
        kkt_error += hdiff * hdiff;
      }
      ++lift_index;
    }
  }
  kkt_error += sto_constraints_->KKTError();
  return kkt_error;
}


inline double SwitchingTimeOptimization::totalCost() const {
  return cost_val_;
}


inline void SwitchingTimeOptimization::integrateSolution(
    const OCP& ocp, std::shared_ptr<ContactSequence>& contact_sequence,
    const double primal_step_size, const double dual_step_size, 
    const Direction& d) const {
  if (is_sto_enabled_) {
    const int N_impulse = ocp.discrete().N_impulse();
    for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
      if (ocp.discrete().isSTOEnabledImpulse(impulse_index)) {
        const double ts = contact_sequence->impulseTime(impulse_index);
        const double ts_new = ts + primal_step_size * d.aux[impulse_index].dts;
        contact_sequence->setImpulseTime(impulse_index, ts_new);
      }
    }
    const int N_lift = ocp.discrete().N_lift();
    for (int lift_index=0; lift_index<N_lift; ++lift_index) {
      if (ocp.discrete().isSTOEnabledLift(lift_index)) {
        const double ts = contact_sequence->liftTime(lift_index);
        const double ts_new = ts + primal_step_size * d.lift[lift_index].dts;
        contact_sequence->setLiftTime(lift_index, ts_new);
      }
    }
    sto_constraints_->updateSlack(primal_step_size);
    sto_constraints_->updateDual(dual_step_size);
  }
}


inline void SwitchingTimeOptimization::setSTORegularization(
    const STORegularization& sto_reg) {
  sto_reg_ = sto_reg;
}


inline void SwitchingTimeOptimization::setSTORegularization(
    const STORegularizationType& reg_type, const double w) {
  sto_reg_ = STORegularization(reg_type, w);
}

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HXX_ 