#include "robotoc/hybrid/switching_time_optimization.hpp"

#include <cassert>


namespace robotoc {

SwitchingTimeOptimization::SwitchingTimeOptimization(const OCP& ocp) 
  : sto_cost_(ocp.sto_cost()), 
    sto_constraints_(ocp.sto_constraints()),
    sto_reg_(0),
    kkt_error_(0),
    cost_val_(0),
    h_phase_(Eigen::VectorXd::Zero(2*ocp.reservedNumDiscreteEvents()+1)),
    reserved_num_switches_(2*ocp.reservedNumDiscreteEvents()),
    is_sto_enabled_(ocp.isSTOEnabled()) {
}


SwitchingTimeOptimization::SwitchingTimeOptimization() 
  : sto_cost_(), 
    sto_constraints_(),
    sto_reg_(0),
    kkt_error_(0),
    cost_val_(0),
    h_phase_(),
    reserved_num_switches_(0),
    is_sto_enabled_(false) {
}


void SwitchingTimeOptimization::reserve(const OCP& ocp) {
  if (reserved_num_switches_ < 2*ocp.reservedNumDiscreteEvents()) {
    h_phase_.resize(2*ocp.reservedNumDiscreteEvents()+1);
    h_phase_.setZero();
    reserved_num_switches_ = 2*ocp.reservedNumDiscreteEvents();
  }
}


void SwitchingTimeOptimization::setRegularization(const double sto_reg) {
  sto_reg_ = sto_reg;
}


void SwitchingTimeOptimization::initConstraints(const OCP& ocp) {
  if (is_sto_enabled_) {
    reserve(ocp);
    sto_constraints_->setSlack(ocp.discrete());
  }
}


void SwitchingTimeOptimization::computeKKTResidual(
    const OCP& ocp, KKTResidual& kkt_residual) {
  const auto& discretization = ocp.discrete();
  if (is_sto_enabled_) {
    reserve(ocp);
    cost_val_ = sto_cost_->linearizeCost(discretization, kkt_residual); 
    sto_constraints_->linearizeConstraints(discretization, kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
  }
}


void SwitchingTimeOptimization::computeKKTSystem(
    const OCP& ocp, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const auto& discretization = ocp.discrete();
  if (is_sto_enabled_) {
    reserve(ocp);
    cost_val_ = sto_cost_->quadratizeCost(discretization, kkt_matrix, 
                                          kkt_residual); 
    sto_constraints_->linearizeConstraints(discretization, kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
    sto_constraints_->condenseSlackAndDual(discretization, kkt_matrix, 
                                           kkt_residual);
  }
}


void SwitchingTimeOptimization::applyRegularization(
    const OCP& ocp, KKTMatrix& kkt_matrix) const {
  const auto& discretization = ocp.discrete();
  for (int i=0; i<discretization.N_impulse(); ++i) {
    kkt_matrix.aux[i].Qtt += sto_reg_;
  }
  for (int i=0; i<discretization.N_lift(); ++i) {
    kkt_matrix.lift[i].Qtt += sto_reg_;
  }
}


void SwitchingTimeOptimization::computeDirection(const OCP& ocp, 
                                                 const Direction& d) {
  if (is_sto_enabled_) {
    sto_constraints_->expandSlackAndDual(ocp.discrete(), d);
  }
}


double SwitchingTimeOptimization::maxPrimalStepSize() const {
  if (is_sto_enabled_) {
    return sto_constraints_->maxPrimalStepSize();
  }
  else {
    return 1.0;
  }
}


double SwitchingTimeOptimization::maxDualStepSize() const {
  if (is_sto_enabled_) {
    return sto_constraints_->maxDualStepSize();
  }
  else {
    return 1.0;
  }
}


double SwitchingTimeOptimization::KKTError() const {
  return kkt_error_;
}


double SwitchingTimeOptimization::KKTError(const OCP& ocp, 
                                           const KKTResidual& kkt_residual) {
  const auto& discretization = ocp.discrete();
  const int N = discretization.N();
  const int N_impulse = discretization.N_impulse();
  const int N_lift = discretization.N_lift();
  const int N_all = N + 1 + 3*N_impulse + N_lift;
  reserve(ocp);
  h_phase_.setZero();
  for (int stage=0; stage<N; ++stage) {
    h_phase_.coeffRef(discretization.contactPhase(stage))
        += kkt_residual[stage].h;
  }
  for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
    h_phase_.coeffRef(discretization.contactPhaseAfterImpulse(impulse_index))
        += kkt_residual.aux[impulse_index].h;
  }
  for (int lift_index=0; lift_index<N_lift; ++lift_index) {
    h_phase_.coeffRef(discretization.contactPhaseAfterLift(lift_index))
        += kkt_residual.lift[lift_index].h;
  }
  double kkt_error = 0;
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<N_impulse+N_lift; ++event_index) {
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        const double hdiff = h_phase_.coeff(event_index) 
                              - h_phase_.coeff(event_index+1);
        kkt_error += hdiff * hdiff;
      }
      ++impulse_index;
    }
    else {
      assert(discretization.eventType(event_index) == DiscreteEventType::Lift);
      if (discretization.isSTOEnabledLift(lift_index)) {
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


double SwitchingTimeOptimization::totalCost() const {
  return cost_val_;
}


void SwitchingTimeOptimization::integrateSolution(OCP& ocp, 
                                                  const double primal_step_size,
                                                  const double dual_step_size, 
                                                  const Direction& d) const {
  auto& contact_sequence = ocp.contact_sequence_nonconst();
  const auto& discretization = ocp.discrete();
  if (is_sto_enabled_) {
    const int N_impulse = discretization.N_impulse();
    for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        const double ts = contact_sequence->impulseTime(impulse_index);
        const double ts_new = ts + primal_step_size * d.aux[impulse_index].dts;
        contact_sequence->setImpulseTime(impulse_index, ts_new);
      }
    }
    const int N_lift = discretization.N_lift();
    for (int lift_index=0; lift_index<N_lift; ++lift_index) {
      if (discretization.isSTOEnabledLift(lift_index)) {
        const double ts = contact_sequence->liftTime(lift_index);
        const double ts_new = ts + primal_step_size * d.lift[lift_index].dts;
        contact_sequence->setLiftTime(lift_index, ts_new);
      }
    }
    sto_constraints_->updateSlack(primal_step_size);
    sto_constraints_->updateDual(dual_step_size);
  }
}

} // namespace robotoc