#ifndef ROBOTOC_OCP_HXX_
#define ROBOTOC_OCP_HXX_

#include <stdexcept>
#include <cassert>

#include "robotoc/ocp/ocp.hpp"


namespace robotoc {

inline OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
                const std::shared_ptr<Constraints>& constraints, 
                const std::shared_ptr<STOCostFunction>& sto_cost, 
                const std::shared_ptr<STOConstraints>& sto_constraints, 
                const double T, const int N, 
                const int max_num_each_discrete_events) 
  : data(N, SplitOCP(robot, cost, constraints)), 
    aux(max_num_each_discrete_events, SplitOCP(robot, cost, constraints)),
    lift(max_num_each_discrete_events, SplitOCP(robot, cost, constraints)),
    impulse(max_num_each_discrete_events, ImpulseSplitOCP(robot, cost, constraints)),
    terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(sto_cost),
    sto_constraints_(sto_constraints),
    discretization_(T, N, max_num_each_discrete_events),
    T_(T),
    N_(N),
    max_num_each_discrete_events_(max_num_each_discrete_events),
    is_sto_enabled_(true) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_each_discrete_events < 0) {
      throw std::out_of_range(
          "invalid argument: max_num_each_discrete_events must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  discretization_.setDiscretizationMethod(DiscretizationMethod::PhaseBased);
}


inline OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
                const std::shared_ptr<Constraints>& constraints, 
                const double T, const int N,
                const int max_num_each_discrete_events) 
  : data(N, SplitOCP(robot, cost, constraints)), 
    aux(max_num_each_discrete_events, SplitOCP(robot, cost, constraints)),
    lift(max_num_each_discrete_events, SplitOCP(robot, cost, constraints)),
    impulse(max_num_each_discrete_events, ImpulseSplitOCP(robot, cost, constraints)),
    terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(),
    sto_constraints_(),
    discretization_(T, N, max_num_each_discrete_events),
    T_(T),
    N_(N),
    max_num_each_discrete_events_(max_num_each_discrete_events),
    is_sto_enabled_(false) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_each_discrete_events < 0) {
      throw std::out_of_range(
          "invalid argument: max_num_each_discrete_events must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline OCP::OCP() 
  : data(), 
    aux(),
    lift(),
    impulse(),
    terminal(),
    robot_(), 
    cost_(),
    constraints_(),
    sto_cost_(),
    sto_constraints_(),
    discretization_(),
    T_(0),
    N_(0),
    is_sto_enabled_(false) {
}


inline void OCP::resize(const int max_num_each_discrete_events) {
  assert(max_num_each_discrete_events >= 0);
  assert(impulse.size() == discretization_.maxNumEachDiscreteEvents());
  assert(aux.size() == discretization_.maxNumEachDiscreteEvents());
  assert(lift.size() == discretization_.maxNumEachDiscreteEvents());
  if (max_num_each_discrete_events > max_num_each_discrete_events_) {
    const auto discretization_method = discretization_.discretizationMethod();
    discretization_ = TimeDiscretization(T_, N_, max_num_each_discrete_events);
    discretization_.setDiscretizationMethod(discretization_method);
    while (max_num_each_discrete_events > impulse.size()) {
      impulse.emplace_back(robot_, cost_, constraints_);
      aux.emplace_back(robot_, cost_, constraints_);
    }
    while (max_num_each_discrete_events > lift.size()) {
      lift.emplace_back(robot_, cost_, constraints_);
    }
    max_num_each_discrete_events_ = max_num_each_discrete_events;
  }
}


inline void OCP::setDiscretizationMethod(
    const DiscretizationMethod discretization_method) {
  if (!is_sto_enabled_) {
    // if is_sto_enabled_ is ture, the discretization method is fixed to 
    // DiscretizationMethod::PhaseBased.
    discretization_.setDiscretizationMethod(discretization_method);
  }
}


inline void OCP::discretize(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  discretization_.discretize(contact_sequence, t);
}


inline void OCP::meshRefinement(
    const std::shared_ptr<ContactSequence>& contact_sequence, const double t) {
  discretization_.meshRefinement(contact_sequence, t);
}


inline const TimeDiscretization& OCP::discrete() const {
  return discretization_;
}


inline const Robot& OCP::robot() const {
  return robot_;
}


inline const std::shared_ptr<CostFunction>& OCP::cost() const {
  return cost_;
}


inline const std::shared_ptr<Constraints>& OCP::constraints() const {
  return constraints_;
}


inline const std::shared_ptr<STOCostFunction>& OCP::sto_cost() const {
  return sto_cost_;
}


inline const std::shared_ptr<STOConstraints>& OCP::sto_constraints() const {
  return sto_constraints_;
}


inline double OCP::T() const {
  return T_;
}


inline int OCP::N() const {
  return N_;
}


inline int OCP::maxNumEachDiscreteEvents() const {
  return max_num_each_discrete_events_;
}


inline bool OCP::isSTOEnabled() const {
  return is_sto_enabled_;
}


inline void OCP::disp(std::ostream& os) const {
  os << "OCP: " << std::endl;
  os << "T: " << T_ << std::endl;
  os << "N: " << N_ << std::endl;
  os << "max_num_each_discrete_events: " << max_num_each_discrete_events_ << std::endl;
  os << robot_ << std::endl;
  os << discretization_ << std::endl;
}


inline std::ostream& operator<<(std::ostream& os, const OCP& ocp) {
  ocp.disp(os);
  return os;
}

} // namespace robotoc

#endif // ROBOTOC_OCP_HXX_ 