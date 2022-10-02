#include "robotoc/ocp/ocp.hpp"

#include <stdexcept>


namespace robotoc {

OCP::OCP(const Robot& _robot, const std::shared_ptr<CostFunction>& _cost, 
         const std::shared_ptr<Constraints>& _constraints, 
         const std::shared_ptr<STOCostFunction>& _sto_cost, 
         const std::shared_ptr<STOConstraints>& _sto_constraints, 
         const std::shared_ptr<ContactSequence>& _contact_sequence,
         const double _T, const int _N, const int _reserved_num_discrete_events) 
  : robot(_robot),
    cost(_cost),
    constraints(_constraints),
    sto_cost(_sto_cost),
    sto_constraints(_sto_constraints),
    contact_sequence(_contact_sequence),
    T(_T),
    N(_N),
    reserved_num_discrete_events(_reserved_num_discrete_events) {
  if (_T <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'T' must be positive!");
  }
  if (_N <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'N' must be positive!");
  }
  if (_reserved_num_discrete_events < 0) {
    throw std::out_of_range("[OCP] invalid argument: 'reserved_num_discrete_events' must be non-negative!");
  }
}


OCP::OCP(const Robot& _robot, const std::shared_ptr<CostFunction>& _cost, 
         const std::shared_ptr<Constraints>& _constraints, 
         const std::shared_ptr<ContactSequence>& _contact_sequence,
         const double _T, const int _N, const int _reserved_num_discrete_events) 
  : robot(_robot),
    cost(_cost),
    constraints(_constraints),
    sto_cost(nullptr),
    sto_constraints(nullptr),
    contact_sequence(_contact_sequence),
    T(_T),
    N(_N),
    reserved_num_discrete_events(_reserved_num_discrete_events) {
  if (_T <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'T' must be positive!");
  }
  if (_N <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'N' must be positive!");
  }
  if (_reserved_num_discrete_events < 0) {
    throw std::out_of_range("[OCP] invalid argument: 'reserved_num_discrete_events' must be non-negative!");
  }
}


OCP::OCP(const Robot& _robot, const std::shared_ptr<CostFunction>& _cost, 
         const std::shared_ptr<Constraints>& _constraints, 
         const double _T, const int _N) 
  : robot(_robot),
    cost(_cost),
    constraints(_constraints),
    sto_cost(nullptr),
    sto_constraints(nullptr),
    contact_sequence(nullptr),
    T(_T),
    N(_N),
    reserved_num_discrete_events(0) {
  if (_T <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'T' must be positive!");
  }
  if (_N <= 0) {
    throw std::out_of_range("[OCP] invalid argument: 'N' must be positive!");
  }
}


OCP::OCP() 
  : robot(), 
    cost(nullptr),
    constraints(nullptr),
    sto_cost(nullptr),
    sto_constraints(nullptr),
    contact_sequence(nullptr),
    T(0),
    N(0),
    reserved_num_discrete_events(0) {
}


void OCP::disp(std::ostream& os) const {
  os << "OCP: " << std::endl;
  os << "  T: " << T << std::endl;
  os << "  N: " << N << std::endl;
  os << "  reserved_num_discrete_events: " << reserved_num_discrete_events << std::endl;
  os << robot << std::endl;
}


std::ostream& operator<<(std::ostream& os, const OCP& ocp) {
  ocp.disp(os);
  return os;
}

} // namespace robotoc