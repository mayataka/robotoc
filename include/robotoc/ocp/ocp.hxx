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
                const std::shared_ptr<ContactSequence>& contact_sequence,
                const double T, const int N) 
  : data(N, SplitOCP(robot, cost, constraints)), 
    aux(contact_sequence->reservedNumDiscreteEvents(), SplitOCP(robot, cost, constraints)),
    lift(contact_sequence->reservedNumDiscreteEvents(), SplitOCP(robot, cost, constraints)),
    impulse(contact_sequence->reservedNumDiscreteEvents(), ImpulseSplitOCP(robot, cost, constraints)),
    terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(sto_cost),
    sto_constraints_(sto_constraints),
    contact_sequence_(contact_sequence),
    discretization_(T, N, contact_sequence->reservedNumDiscreteEvents()),
    T_(T),
    N_(N),
    reserved_num_discrete_events_(contact_sequence->reservedNumDiscreteEvents()),
    is_sto_enabled_(true) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    const auto& min_dt = sto_constraints->minimumDwellTimes();
    double sum_min_dt = 0;
    for (const auto e : min_dt) {
      sum_min_dt += e;
    }
    if (T <= sum_min_dt) {
      throw std::out_of_range(
          "invalid argument: sum of the minimum dwell-times must be smaller than T!");
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
                const std::shared_ptr<ContactSequence>& contact_sequence,
                const double T, const int N) 
  : data(N, SplitOCP(robot, cost, constraints)), 
    aux(contact_sequence->reservedNumDiscreteEvents(), SplitOCP(robot, cost, constraints)),
    lift(contact_sequence->reservedNumDiscreteEvents(), SplitOCP(robot, cost, constraints)),
    impulse(contact_sequence->reservedNumDiscreteEvents(), ImpulseSplitOCP(robot, cost, constraints)),
    terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(),
    sto_constraints_(),
    contact_sequence_(contact_sequence),
    discretization_(T, N, contact_sequence->reservedNumDiscreteEvents()),
    T_(T),
    N_(N),
    reserved_num_discrete_events_(contact_sequence->reservedNumDiscreteEvents()),
    is_sto_enabled_(false) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
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
    contact_sequence_(),
    discretization_(),
    T_(0),
    N_(0),
    reserved_num_discrete_events_(0),
    is_sto_enabled_(false) {
}


inline void OCP::reserve() {
  assert(impulse.size() == reserved_num_discrete_events_);
  assert(aux.size() == reserved_num_discrete_events_);
  assert(lift.size() == reserved_num_discrete_events_);
  const int new_reserved_num_discrete_events 
      = discretization_.reservedNumDiscreteEvents();
  if (new_reserved_num_discrete_events > reserved_num_discrete_events_) {
    while (impulse.size() < new_reserved_num_discrete_events) {
      impulse.emplace_back(robot_, cost_, constraints_);
    }
    while (aux.size() < new_reserved_num_discrete_events) {
      aux.emplace_back(robot_, cost_, constraints_);
    }
    while (lift.size() < new_reserved_num_discrete_events) {
      lift.emplace_back(robot_, cost_, constraints_);
    }
    reserved_num_discrete_events_ = new_reserved_num_discrete_events;
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


inline void OCP::discretize(const double t) {
  discretization_.discretize(contact_sequence_, t);
  reserve();
}


inline void OCP::meshRefinement(const double t) {
  discretization_.meshRefinement(contact_sequence_, t);
  reserve();
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


inline const std::shared_ptr<ContactSequence>& OCP::contact_sequence() const {
  return contact_sequence_;
}


inline double OCP::T() const {
  return T_;
}


inline int OCP::N() const {
  return N_;
}


inline int OCP::reservedNumDiscreteEvents() const {
  return reserved_num_discrete_events_;
}


inline bool OCP::isSTOEnabled() const {
  return is_sto_enabled_;
}


inline void OCP::disp(std::ostream& os) const {
  os << "OCP: " << std::endl;
  os << "T: " << T_ << std::endl;
  os << "N: " << N_ << std::endl;
  os << "reserved_num_discrete_events: " << reserved_num_discrete_events_ << std::endl;
  os << robot_ << std::endl;
  os << discretization_ << std::endl;
}


inline std::ostream& operator<<(std::ostream& os, const OCP& ocp) {
  ocp.disp(os);
  return os;
}

} // namespace robotoc

#endif // ROBOTOC_OCP_HXX_ 