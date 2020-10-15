#include "idocp/cost/impulse_force_cost.hpp"

#include <iostream>


namespace idocp {

ImpulseForceCost::ImpulseForceCost(const Robot& robot)
  : ImpulseCostFunctionComponentBase(),
    max_point_contacts_(robot.max_point_contacts()),
    max_dimf_(robot.max_dimf()) {
}


ImpulseForceCost::ImpulseForceCost()
  : ImpulseCostFunctionComponentBase(),
    max_point_contacts_(0),
    max_dimf_(0) {
}


ImpulseForceCost::~ImpulseForceCost() {
}



void ImpulseForceCost::set_f_ref(const std::vector<Eigen::Vector3d>& f_ref) {
  if (f_ref.size() == max_point_contacts_) {
    f_ref_ = f_ref;
  }
  else {
    std::cout << "invalid argment in set_f_ref(): size of f_ref must be " 
              << max_point_contacts_ << std::endl;
  }
}


void ImpulseForceCost::set_f_weight(
    const std::vector<Eigen::Vector3d>& f_weight) {
  if (f_weight.size() == max_point_contacts_) {
    f_weight_ = f_weight;
  }
  else {
    std::cout << "invalid argment in set_f_weight(): size of f_ref must be " 
              << max_point_contacts_ << std::endl;
  }
}


double ImpulseForceCost::l(Robot& robot, CostFunctionData& data, const double t, 
                           const ImpulseSplitSolution& s) const {
  double l = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      l += (f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array()) 
                                 * (s.f[i].array()-f_ref_[i].array())).sum();
    }
  }
  return 0.5 * l;
}


void ImpulseForceCost::lf(Robot& robot, CostFunctionData& data, const double t, 
                          const ImpulseSplitSolution& s, 
                          ImpulseKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array());
      dimf_stack += 3;
    }
  }
}


void ImpulseForceCost::lff(Robot& robot, CostFunctionData& data, const double t, 
                           const ImpulseSplitSolution& s, 
                           ImpulseKKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
          += f_weight_[i];
      dimf_stack += 3;
    }
  }
}

} // namespace idocp