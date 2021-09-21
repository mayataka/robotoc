#ifndef IDOCP_SWITCHING_TIME_CONSTRAINTS_HXX_
#define IDOCP_SWITCHING_TIME_CONSTRAINTS_HXX_

#include "idocp/hybrid/switching_time_constraints.hpp"


namespace idocp {

inline SwitchingTimeConstraints::SwitchingTimeConstraints(
    const int max_num_switches, const double min_dt, const double min_dt0, 
    const double min_dtf, const double barrier, 
    const double fraction_to_boundary_rule) 
  : dtlb_(max_num_switches+1, DwellTimeLowerBound(barrier, 
                                                  fraction_to_boundary_rule)),
    min_dt_(min_dt), 
    min_dt0_(min_dt0), 
    min_dtf_(min_dtf),
    num_switches_(0),
    primal_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)), 
    dual_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)) {
}


inline SwitchingTimeConstraints::SwitchingTimeConstraints()
  : dtlb_(),
    min_dt_(0), 
    min_dt0_(0), 
    min_dtf_(0),
    num_switches_(0),
    primal_step_size_(), 
    dual_step_size_() {
}


inline SwitchingTimeConstraints::~SwitchingTimeConstraints() {
}


inline void SwitchingTimeConstraints::setSlack(
    const HybridOCPDiscretization& discretization) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    dtlb_[0].setSlack(min_dt0_, discretization.t(0), discretization.t_impulse(0));
  }
  else {
    dtlb_[0].setSlack(min_dt0_, discretization.t(0), discretization.t_lift(0));
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_,
                                      discretization.t_impulse(impulse_index),  
                                      discretization.t_impulse(impulse_index+1));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_, 
                                      discretization.t_impulse(impulse_index),  
                                      discretization.t_lift(lift_index));
      }
      ++impulse_index;
    }
    else {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_,
                                      discretization.t_lift(lift_index),  
                                      discretization.t_impulse(impulse_index));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_,
                                      discretization.t_lift(lift_index),  
                                      discretization.t_lift(lift_index+1));
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    dtlb_[num_events].setSlack(min_dtf_, discretization.t_impulse(impulse_index), 
                               discretization.t(discretization.N()));
  }
  else {
    dtlb_[num_events].setSlack(min_dtf_, discretization.t_lift(lift_index), 
                               discretization.t(discretization.N()));
  }
}


inline void SwitchingTimeConstraints::evalConstraint(
    const HybridOCPDiscretization& discretization) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_impulse(0));
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_lift(0));
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_impulse(impulse_index+1));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_lift(lift_index));
        }
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_impulse(impulse_index));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_lift(lift_index+1));
        }
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_impulse(impulse_index), 
                                       discretization.t(discretization.N()));
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_lift(lift_index), 
                                       discretization.t(discretization.N()));
    }
  }
}


inline void SwitchingTimeConstraints::linearizeConstraints(
    const HybridOCPDiscretization& discretization, 
    KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_impulse(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.aux[0]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_lift(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.lift[0]);
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_impulse(impulse_index+1));
          dtlb_[event_index+1].evalDerivatives_lub( 
              kkt_residual.aux[impulse_index], 
              kkt_residual.aux[impulse_index+1]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_lift(lift_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.aux[impulse_index], kkt_residual.lift[lift_index]);
        }
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_impulse(impulse_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.aux[impulse_index]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_lift(lift_index+1));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.lift[lift_index+1]);
        }
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_impulse(impulse_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_lift(lift_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.lift[lift_index]);
    }
  }
}


inline void SwitchingTimeConstraints::condenseSlackAndDual(
    const HybridOCPDiscretization& discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_impulse(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.aux[0]);
      dtlb_[0].condenseSlackAndDual_lb(kkt_matrix.aux[0], kkt_residual.aux[0]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt0_, discretization.t(0), 
                              discretization.t_lift(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.lift[0]);
      dtlb_[0].condenseSlackAndDual_lb(kkt_matrix.lift[0], kkt_residual.lift[0]);
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_impulse(impulse_index+1));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.aux[impulse_index], kkt_residual.aux[impulse_index+1]);
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
              kkt_matrix.aux[impulse_index+1], kkt_residual.aux[impulse_index+1]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_impulse(impulse_index), 
              discretization.t_lift(lift_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.aux[impulse_index], kkt_residual.lift[lift_index]);
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
        }
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_impulse(impulse_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.aux[impulse_index]);
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_, discretization.t_lift(lift_index), 
              discretization.t_lift(lift_index+1));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.lift[lift_index+1]);
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index],
              kkt_matrix.lift[lift_index+1], kkt_residual.lift[lift_index+1]);
        }
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_impulse(impulse_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.aux[impulse_index]);
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dtf_, 
                                       discretization.t_lift(lift_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.lift[lift_index]);
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


inline void SwitchingTimeConstraints::expandSlackAndDual(
    const HybridOCPDiscretization& discretization, const Direction& d) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.aux[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
    else {
      primal_step_size_.coeffRef(0) = 1.0;
      dual_step_size_.coeffRef(0) = 1.0;
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.lift[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
    else {
      primal_step_size_.coeffRef(0) = 1.0;
      dual_step_size_.coeffRef(0) = 1.0;
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].expandSlackAndDual_lub((-d.aux[impulse_index].dts),
                                                      (-d.aux[impulse_index+1].dts));
        }
        else {
          dtlb_[event_index+1].expandSlackAndDual_lub((-d.aux[impulse_index].dts),
                                                      (-d.lift[lift_index].dts));
        }
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      else {
        primal_step_size_.coeffRef(event_index+1) = 1.0;
        dual_step_size_.coeffRef(event_index+1) = 1.0;
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].expandSlackAndDual_lub((-d.lift[lift_index].dts),
                                                      (-d.aux[impulse_index+1].dts));
        }
        else {
          dtlb_[event_index+1].expandSlackAndDual_lub((-d.lift[lift_index].dts),
                                                      (-d.lift[lift_index+1].dts));
        }
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      else {
        primal_step_size_.coeffRef(event_index+1) = 1.0;
        dual_step_size_.coeffRef(event_index+1) = 1.0;
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(-d.aux[impulse_index].dts);
      primal_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxDualStepSize();
    }
    else {
      primal_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxDualStepSize();
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(-d.lift[lift_index].dts);
      primal_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxDualStepSize();
    }
    else {
      primal_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) 
          = dtlb_[num_events].maxDualStepSize();
    }
  }
}


inline double SwitchingTimeConstraints::maxPrimalStepSize() const {
  if (num_switches_ > 0) {
    return primal_step_size_.head(num_switches_).minCoeff();
  }
  else {
    return 1.0;
  }
}


inline double SwitchingTimeConstraints::maxDualStepSize() const {
  if (num_switches_ > 0) {
    return dual_step_size_.head(num_switches_).minCoeff();
  }
  else {
    return 1.0;
  }
}


inline void SwitchingTimeConstraints::updateSlack(const double step_size) {
  for (int i=0; i<num_switches_; ++i) {
    dtlb_[i].updateSlack(step_size);
  }
}


inline void SwitchingTimeConstraints::updateDual(const double step_size) {
  for (int i=0; i<num_switches_; ++i) {
    dtlb_[i].updateDual(step_size);
  }
}

} // namespace idocp

#endif // IDOCP_SWITCHING_TIME_CONSTRAINTS_HXX_ 