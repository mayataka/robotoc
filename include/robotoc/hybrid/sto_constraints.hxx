#ifndef ROBOTOC_STO_CONSTRAINTS_HXX_ 
#define ROBOTOC_STO_CONSTRAINTS_HXX_

#include "robotoc/hybrid/sto_constraints.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

inline STOConstraints::STOConstraints(const int max_num_switches, 
                                      const double min_dt, const double barrier, 
                                      const double fraction_to_boundary_rule) 
  : dtlb_(max_num_switches+1, DwellTimeLowerBound(barrier, 
                                                  fraction_to_boundary_rule)),
    min_dt_(max_num_switches+1, min_dt), 
    max_num_switches_(max_num_switches),
    num_switches_(0),
    primal_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)), 
    dual_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)) {
  try {
    if (min_dt < 0) {
      throw std::out_of_range(
          "Invalid argment: min_dt must be non-negative!");
    }
    if (barrier <= 0) {
      throw std::out_of_range(
          "Invalid argment: barrirer must be positive!");
    }
    if (fraction_to_boundary_rule <= 0) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be positive!");
    }
    if (fraction_to_boundary_rule >= 1) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be less than 1!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline STOConstraints::STOConstraints(const int max_num_switches, 
                                      const std::vector<double>& min_dt, 
                                      const double barrier, 
                                      const double fraction_to_boundary_rule) 
  : dtlb_(max_num_switches+1, DwellTimeLowerBound(barrier, 
                                                  fraction_to_boundary_rule)),
    min_dt_(min_dt), 
    max_num_switches_(max_num_switches),
    num_switches_(0),
    primal_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)), 
    dual_step_size_(Eigen::VectorXd::Zero(max_num_switches+1)) {
  try {
    for (const auto e : min_dt) {
      if (e < 0.) {
        throw std::out_of_range(
            "Invalid argment: min_dt must be non-negative!");
      }
    }
    if (barrier <= 0) {
      throw std::out_of_range(
          "Invalid argment: barrirer must be positive!");
    }
    if (fraction_to_boundary_rule <= 0) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be positive!");
    }
    if (fraction_to_boundary_rule >= 1) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be less than 1!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  while (min_dt_.size() < (max_num_switches+1)) {
    min_dt_.push_back(k_min_dt);
  }
}


inline STOConstraints::STOConstraints()
  : dtlb_(),
    min_dt_(0), 
    num_switches_(0),
    primal_step_size_(), 
    dual_step_size_() {
}


inline STOConstraints::~STOConstraints() {
}


inline void STOConstraints::setSlack(
    const HybridOCPDiscretization& discretization) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    dtlb_[0].setSlack(min_dt_[0], discretization.t(0), discretization.t_impulse(0));
  }
  else {
    dtlb_[0].setSlack(min_dt_[0], discretization.t(0), discretization.t_lift(0));
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      discretization.t_impulse(impulse_index),  
                                      discretization.t_impulse(impulse_index+1));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1], 
                                      discretization.t_impulse(impulse_index),  
                                      discretization.t_lift(lift_index));
      }
      ++impulse_index;
    }
    else {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      discretization.t_lift(lift_index),  
                                      discretization.t_impulse(impulse_index));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      discretization.t_lift(lift_index),  
                                      discretization.t_lift(lift_index+1));
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    dtlb_[num_events].setSlack(min_dt_[num_events], 
                               discretization.t_impulse(impulse_index), 
                               discretization.t(discretization.N()));
  }
  else {
    dtlb_[num_events].setSlack(min_dt_[num_events], 
                               discretization.t_lift(lift_index), 
                               discretization.t(discretization.N()));
  }
}


inline void STOConstraints::evalConstraint(
    const HybridOCPDiscretization& discretization) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], discretization.t(0), 
                              discretization.t_impulse(0));
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], discretization.t(0), 
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
              min_dt_[event_index+1], discretization.t_impulse(impulse_index), 
              discretization.t_impulse(impulse_index+1));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], discretization.t_impulse(impulse_index), 
              discretization.t_lift(lift_index));
        }
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], discretization.t_lift(lift_index), 
              discretization.t_impulse(impulse_index));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], discretization.t_lift(lift_index), 
              discretization.t_lift(lift_index+1));
        }
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       discretization.t_impulse(impulse_index), 
                                       discretization.t(discretization.N()));
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       discretization.t_lift(lift_index), 
                                       discretization.t(discretization.N()));
    }
  }
}


inline void STOConstraints::linearizeConstraints(
    const HybridOCPDiscretization& discretization, 
    KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], discretization.t(0), 
                              discretization.t_impulse(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.aux[0]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], discretization.t(0), 
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
              min_dt_[event_index+1], discretization.t_impulse(impulse_index), 
              discretization.t_impulse(impulse_index+1));
          dtlb_[event_index+1].evalDerivatives_lub( 
              kkt_residual.aux[impulse_index], 
              kkt_residual.aux[impulse_index+1]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], discretization.t_impulse(impulse_index), 
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
               min_dt_[event_index+1], discretization.t_lift(lift_index), 
              discretization.t_impulse(impulse_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.aux[impulse_index]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], discretization.t_lift(lift_index), 
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
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       discretization.t_impulse(impulse_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       discretization.t_lift(lift_index), 
                                       discretization.t(discretization.N()));
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.lift[lift_index]);
    }
  }
}


inline void STOConstraints::condenseSlackAndDual(
    const HybridOCPDiscretization& discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].condenseSlackAndDual_lb(kkt_matrix.aux[0], kkt_residual.aux[0]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
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
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
              kkt_matrix.aux[impulse_index+1], kkt_residual.aux[impulse_index+1]);
        }
        else {
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
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
        }
        else {
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
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


inline void STOConstraints::expandSlackAndDual(
    const HybridOCPDiscretization& discretization, const Direction& d) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  primal_step_size_.fill(1.0);
  dual_step_size_.fill(1.0);
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.aux[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
  }
  else {
    if (discretization.isSTOEnabledLift(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.lift[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (discretization.isSTOEnabledImpulse(impulse_index)) {
        dtlb_[event_index+1].expandSlackAndDual_lub(d.aux[impulse_index].dts,
                                                    d.aux[impulse_index].dts_next);
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      ++impulse_index;
    }
    else {
      if (discretization.isSTOEnabledLift(lift_index)) {
        dtlb_[event_index+1].expandSlackAndDual_lub(d.lift[lift_index].dts,
                                                    d.lift[lift_index].dts_next);
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      ++lift_index;
    }
  }
  if (discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(d.aux[impulse_index].dts);
      primal_step_size_.coeffRef(num_events) = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) = dtlb_[num_events].maxDualStepSize();
    } 
  }
  else {
    if (discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(d.lift[lift_index].dts);
      primal_step_size_.coeffRef(num_events) = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) = dtlb_[num_events].maxDualStepSize();
    }
  }
}


inline double STOConstraints::maxPrimalStepSize() const {
  if (num_switches_ > 0) {
    return primal_step_size_.head(num_switches_+1).minCoeff();
  }
  else {
    return 1.0;
  }
}


inline double STOConstraints::maxDualStepSize() const {
  if (num_switches_ > 0) {
    return dual_step_size_.head(num_switches_+1).minCoeff();
  }
  else {
    return 1.0;
  }
}


inline void STOConstraints::updateSlack(const double step_size) {
  for (int i=0; i<num_switches_+1; ++i) {
    dtlb_[i].updateSlack(step_size);
  }
}


inline void STOConstraints::updateDual(const double step_size) {
  for (int i=0; i<num_switches_+1; ++i) {
    dtlb_[i].updateDual(step_size);
  }
}


inline double STOConstraints::KKTError() const {
  double err = 0;
  for (int i=0; i<num_switches_+1; ++i) {
    err += dtlb_[i].KKTError();
  }
  return err;
}


inline void STOConstraints::setBarrier(const double barrier) {
  for (auto& e : dtlb_) {
    e.setBarrier(barrier);
  }
}


inline void STOConstraints::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  for (auto& e : dtlb_) {
    e.setFractionToBoundaryRule(fraction_to_boundary_rule);
  }
}


inline void STOConstraints::setMinimumDwellTimes(const double min_dt) {
  try {
    if (min_dt < 0) {
      throw std::out_of_range(
          "Invalid argment: min_dt must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  for (auto& e : min_dt_) {
    e = min_dt;
  }
}


inline void STOConstraints::setMinimumDwellTimes(
    const std::vector<double>& min_dt) {
  min_dt_ = min_dt;
  while (min_dt_.size() < (max_num_switches_+1)) {
    min_dt_.push_back(k_min_dt);
  }
}

} // namespace robotoc

#endif // ROBOTOC_STO_CONSTRAINTS_HXX_ 