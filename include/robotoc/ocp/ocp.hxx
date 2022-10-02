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
  : terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(sto_cost),
    sto_constraints_(sto_constraints),
    contact_sequence_(contact_sequence),
    T_(T),
    N_(N),
    reserved_num_discrete_events_(contact_sequence->reservedNumDiscreteEvents()),
    is_sto_enabled_(true) {
  if (T <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'T' must be positive!");
  }
  if (N <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'N' must be positive!");
  }
  const auto& min_dt = sto_constraints->getMinimumDwellTimes();
  double sum_min_dt = 0;
  for (const auto e : min_dt) {
    sum_min_dt += e;
  }
  if (T <= sum_min_dt) {
    throw std::out_of_range(
        "[OCP] invalid argument: sum of the minimum dwell-times must be smaller than T!");
  }
}


inline OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
                const std::shared_ptr<Constraints>& constraints, 
                const std::shared_ptr<ContactSequence>& contact_sequence,
                const double T, const int N) 
  : terminal(TerminalOCP(robot, cost, constraints)),
    robot_(robot),
    cost_(cost),
    constraints_(constraints),
    sto_cost_(),
    sto_constraints_(),
    contact_sequence_(contact_sequence),
    T_(T),
    N_(N),
    reserved_num_discrete_events_(contact_sequence->reservedNumDiscreteEvents()),
    is_sto_enabled_(false) {
  if (T <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'T' must be positive!");
  }
  if (N <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'N' must be positive!");
  }
}


inline OCP::OCP() 
  : terminal(),
    robot_(), 
    cost_(),
    constraints_(),
    sto_cost_(),
    sto_constraints_(),
    contact_sequence_(),
    T_(0),
    N_(0),
    reserved_num_discrete_events_(0),
    is_sto_enabled_(false) {
}


inline void OCP::reserve() {
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
}


inline std::ostream& operator<<(std::ostream& os, const OCP& ocp) {
  ocp.disp(os);
  return os;
}

} // namespace robotoc

#endif // ROBOTOC_OCP_HXX_ 