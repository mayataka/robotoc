#include "idocp/cost/impulse_force_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

ImpulseForceCost::ImpulseForceCost(const Robot& robot)
  : ImpulseCostFunctionComponentBase(),
    max_point_contacts_(robot.maxPointContacts()),
    max_dimf_(robot.max_dimf()),
    f_ref_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
    f_weight_(robot.maxPointContacts(), Eigen::Vector3d::Zero()) {
}


ImpulseForceCost::ImpulseForceCost()
  : ImpulseCostFunctionComponentBase(),
    max_point_contacts_(0),
    max_dimf_(0),
    f_ref_(),
    f_weight_() {
}


ImpulseForceCost::~ImpulseForceCost() {
}



void ImpulseForceCost::set_f_ref(const std::vector<Eigen::Vector3d>& f_ref) {
  try {
    if (f_ref.size() != max_point_contacts_) {
      throw std::invalid_argument(
          "invalid size: f_ref.size() must be " 
          + std::to_string(max_point_contacts_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  f_ref_ = f_ref;
}


void ImpulseForceCost::set_f_weight(
    const std::vector<Eigen::Vector3d>& f_weight) {
  try {
    if (f_weight.size() != max_point_contacts_) {
      throw std::invalid_argument(
          "invalid size: f_weight.size() must be " 
          + std::to_string(max_point_contacts_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  f_weight_ = f_weight;
}


double ImpulseForceCost::l(Robot& robot, CostFunctionData& data, const double t, 
                           const ImpulseSplitSolution& s) const {
  double l = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isImpulseActive(i)) {
      l += (f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array()) 
                                 * (s.f[i].array()-f_ref_[i].array())).sum();
    }
  }
  return 0.5 * l;
}


void ImpulseForceCost::lf(Robot& robot, CostFunctionData& data, const double t, 
                          const ImpulseSplitSolution& s, 
                          ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isImpulseActive(i)) {
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array());
      dimf_stack += 3;
    }
  }
}


void ImpulseForceCost::lff(Robot& robot, CostFunctionData& data, const double t, 
                           const ImpulseSplitSolution& s, 
                           ImpulseSplitKKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isImpulseActive(i)) {
      kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
          += f_weight_[i];
      dimf_stack += 3;
    }
  }
}

} // namespace idocp