#include "robotoc/sto/switching_time_optimization.hpp"

#include <cassert>


namespace robotoc {

SwitchingTimeOptimization::SwitchingTimeOptimization(const OCP& ocp) 
  : sto_cost_(ocp.sto_cost()), 
    sto_constraints_(ocp.sto_constraints()),
    contact_sequence_(ocp.contact_sequence()),
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
    sto_constraints_->setSlack(ocp.timeDiscretization());
  }
}


void SwitchingTimeOptimization::computeKKTResidual(
    const OCP& ocp, KKTResidual& kkt_residual) {
  if (is_sto_enabled_) {
    reserve(ocp);
    cost_val_ = sto_cost_->linearizeCost(ocp.timeDiscretization(), kkt_residual); 
    sto_constraints_->linearizeConstraints(ocp.timeDiscretization(), kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
  }
}


void SwitchingTimeOptimization::computeKKTSystem(
    const OCP& ocp, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  if (is_sto_enabled_) {
    reserve(ocp);
    cost_val_ = sto_cost_->quadratizeCost(ocp.timeDiscretization(), kkt_matrix, 
                                          kkt_residual); 
    sto_constraints_->linearizeConstraints(ocp.timeDiscretization(), kkt_residual);
    kkt_error_ = KKTError(ocp, kkt_residual);
    sto_constraints_->condenseSlackAndDual(ocp.timeDiscretization(), kkt_matrix, 
                                           kkt_residual);
  }
}


void SwitchingTimeOptimization::applyRegularization(
    const OCP& ocp, KKTMatrix& kkt_matrix) const {
  for (int i=0; i<ocp.timeDiscretization().N_impulse(); ++i) {
    kkt_matrix.aux[i].Qtt += sto_reg_;
  }
  for (int i=0; i<ocp.timeDiscretization().N_lift(); ++i) {
    kkt_matrix.lift[i].Qtt += sto_reg_;
  }
}


void SwitchingTimeOptimization::computeDirection(const OCP& ocp, 
                                                 const Direction& d) {
  if (is_sto_enabled_) {
    sto_constraints_->expandSlackAndDual(ocp.timeDiscretization(), d);
  }
}


// double SwitchingTimeOptimization::maxPrimalStepSize() const {
//   if (is_sto_enabled_) {
//     return sto_constraints_->maxPrimalStepSize();
//   }
//   else {
//     return 1.0;
//   }
// }


// double SwitchingTimeOptimization::maxDualStepSize() const {
//   if (is_sto_enabled_) {
//     return sto_constraints_->maxDualStepSize();
//   }
//   else {
//     return 1.0;
//   }
// }


double SwitchingTimeOptimization::KKTError() const {
  return kkt_error_;
}


double SwitchingTimeOptimization::KKTError(const OCP& ocp, 
                                           const KKTResidual& kkt_residual) {
  const int N = ocp.timeDiscretization().N();
  const int N_impulse = ocp.timeDiscretization().N_impulse();
  const int N_lift = ocp.timeDiscretization().N_lift();
  const int N_all = N + 1 + 3*N_impulse + N_lift;
  reserve(ocp);
  h_phase_.setZero();
  for (int stage=0; stage<N; ++stage) {
    h_phase_.coeffRef(ocp.timeDiscretization().contactPhase(stage))
        += kkt_residual[stage].h;
  }
  for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
    h_phase_.coeffRef(ocp.timeDiscretization().contactPhaseAfterImpulse(impulse_index))
        += kkt_residual.aux[impulse_index].h;
  }
  for (int lift_index=0; lift_index<N_lift; ++lift_index) {
    h_phase_.coeffRef(ocp.timeDiscretization().contactPhaseAfterLift(lift_index))
        += kkt_residual.lift[lift_index].h;
  }
  double kkt_error = 0;
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<N_impulse+N_lift; ++event_index) {
    if (ocp.timeDiscretization().eventType(event_index) == DiscreteEventType::Impulse) {
      if (ocp.timeDiscretization().isSTOEnabledImpulse(impulse_index)) {
        const double hdiff = h_phase_.coeff(event_index) 
                              - h_phase_.coeff(event_index+1);
        kkt_error += hdiff * hdiff;
      }
      ++impulse_index;
    }
    else {
      assert(ocp.timeDiscretization().eventType(event_index) == DiscreteEventType::Lift);
      if (ocp.timeDiscretization().isSTOEnabledLift(lift_index)) {
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


void SwitchingTimeOptimization::integrateSolution(
    const OCP& ocp, std::shared_ptr<ContactSequence>& contact_sequence,
    const double primal_step_size, const double dual_step_size, 
    const Direction& d) const {
  if (is_sto_enabled_) {
    const int N_impulse = ocp.timeDiscretization().N_impulse();
    for (int impulse_index=0; impulse_index<N_impulse; ++impulse_index) {
      if (ocp.timeDiscretization().isSTOEnabledImpulse(impulse_index)) {
        const double ts = contact_sequence->impulseTime(impulse_index);
        const double ts_new = ts + primal_step_size * d.aux[impulse_index].dts;
        contact_sequence->setImpulseTime(impulse_index, ts_new);
      }
    }
    const int N_lift = ocp.timeDiscretization().N_lift();
    for (int lift_index=0; lift_index<N_lift; ++lift_index) {
      if (ocp.timeDiscretization().isSTOEnabledLift(lift_index)) {
        const double ts = contact_sequence->liftTime(lift_index);
        const double ts_new = ts + primal_step_size * d.lift[lift_index].dts;
        contact_sequence->setLiftTime(lift_index, ts_new);
      }
    }
    sto_constraints_->updateSlack(primal_step_size);
    sto_constraints_->updateDual(dual_step_size);
  }
}




void SwitchingTimeOptimization::initConstraints(
    const TimeDiscretization& time_discretization) {
  if (!is_sto_enabled_) return;

  constraint_data_ = sto_constraints_->createConstraintsData(time_discretization);
  sto_constraints_->setSlackAndDual(time_discretization, constraint_data_);
}


bool SwitchingTimeOptimization::isFeasible(
    const TimeDiscretization& time_discretization) {
  if (!is_sto_enabled_) return true;

  return sto_constraints_->isFeasible(time_discretization, constraint_data_);
}


void SwitchingTimeOptimization::evalSTO(
    const TimeDiscretization& time_discretization) {
  if (!is_sto_enabled_) return;

  sto_constraints_->evalConstraint(time_discretization, constraint_data_);
  performance_index_.cost_barrier = constraint_data_.log_barrier;
  performance_index_.primal_feasibility = constraint_data_.primalFeasibility();
}


void SwitchingTimeOptimization::evalKKT(
    const TimeDiscretization& time_discretization,
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  if (!is_sto_enabled_) return;

  const int N = time_discretization.N_grids();
  const int num_discrete_events = time_discretization.grid(N).contact_phase 
                                    - time_discretization.grid(0).contact_phase;
  performance_index_.setZero();
  lt_.setZero(num_discrete_events);
  Qtt_.setZero(num_discrete_events, num_discrete_events);
  Qtt_.diagonal().fill(sto_reg_);
  sto_constraints_->linearizeConstraints(time_discretization, 
                                         constraint_data_, lt_);
  performance_index_.cost_barrier = constraint_data_.log_barrier;
  performance_index_.primal_feasibility = constraint_data_.primalFeasibility();
  performance_index_.dual_feasibility = constraint_data_.dualFeasibility();
  performance_index_.kkt_error = constraint_data_.KKTError();
  performance_index_.cost_barrier = constraint_data_.log_barrier;
  sto_constraints_->condenseSlackAndDual(constraint_data_, lt_, Qtt_);

  int event_index = 0;
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    if (grid.type == GridType::Impulse) {
      kkt_residual[i+1].h -= lt_.coeff(event_index);
      kkt_matrix[i+1].Qtt += Qtt_.coeff(event_index, event_index);
      ++event_index;
    }
    else if (grid.type == GridType::Lift) {
      kkt_residual[i].h -= lt_.coeff(event_index);
      kkt_matrix[i].Qtt += Qtt_.coeff(event_index, event_index);
      ++event_index;
    }
  }

  h_.setZero(num_discrete_events+1);
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    h_.coeffRef(grid.contact_phase) += kkt_residual[i].h;
  }

  event_index = 0;
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    if ((grid.type == GridType::Impulse && time_discretization.grid(i+1).sto)
        || (grid.type == GridType::Lift && grid.sto)) {
      const double hdiff = h_.coeff(event_index) - h_.coeff(event_index+1);
      performance_index_.kkt_error = hdiff * hdiff;
      ++event_index;
    }
  }
}


void SwitchingTimeOptimization::computeStepSizes(
    const TimeDiscretization& time_discretization, Direction& d) {
  if (!is_sto_enabled_) return;

  const int N = time_discretization.N_grids();
  const int num_discrete_events = time_discretization.grid(N).contact_phase 
                                    - time_discretization.grid(0).contact_phase;
  dts_.resize(num_discrete_events);
  int event_index = 0;
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    if (grid.type == GridType::Impulse) {
      dts_.coeffRef(event_index) = d[i+1].dts;
      ++event_index;
    }
    else if (grid.type == GridType::Lift) {
      dts_.coeffRef(event_index) = d[i].dts;
      ++event_index;
    }
  }

  sto_constraints_->expandSlackAndDual(constraint_data_, dts_);
}


double SwitchingTimeOptimization::maxPrimalStepSize() const {
  if (dts_.size() > 0) {
    return sto_constraints_->maxPrimalStepSize();
  }
  else {
    return 1.0;
  }
}


double SwitchingTimeOptimization::maxDualStepSize() const {
  if (dts_.size() > 0) {
    return sto_constraints_->maxDualStepSize();
  }
  else {
    return 1.0;
  }
}


void SwitchingTimeOptimization::integrateSolution(
    const TimeDiscretization& time_discretization, 
    const double primal_step_size, const double dual_step_size, 
    const Direction& d) {
  if (!is_sto_enabled_) return;

  const int N = time_discretization.N_grids();
  const int num_discrete_events = time_discretization.grid(N).contact_phase 
                                    - time_discretization.grid(0).contact_phase;
  for (int i=0; i<N; ++i) {
    const auto& grid = time_discretization.grid(i);
    if (grid.type == GridType::Impulse) {
      const int impulse_index = grid.impulse_index;
      const double ts = contact_sequence_->impulseTime(impulse_index) 
                          + primal_step_size * d[i+1].dts;
      contact_sequence_->setImpulseTime(impulse_index, ts);
    }
    else if (grid.type == GridType::Lift) {
      const int lift_index = grid.lift_index;
      const double ts = contact_sequence_->impulseTime(lift_index)  
                          + primal_step_size * d[i].dts;
      contact_sequence_->setLiftTime(lift_index, ts);
    }
  }
  sto_constraints_->updateSlack(constraint_data_, primal_step_size);
  sto_constraints_->updateDual(constraint_data_, dual_step_size);
}


const PerformanceIndex& SwitchingTimeOptimization::getEval() const {
  return performance_index_;
}

} // namespace robotoc