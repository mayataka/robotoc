#include "idocp/cost/impulse_force_cost.hpp"

#include <iostream>


namespace idocp {

ImpulseForceCost::ImpulseForceCost(const Robot& robot)
  : ImpulseCostFunctionComponentBase(),
    max_dimf_(robot.max_dimf()) {
}


ImpulseForceCost::ImpulseForceCost()
  : ImpulseCostFunctionComponentBase(),
    max_dimf_(0) {
}


ImpulseForceCost::~ImpulseForceCost() {
}



void ImpulseForceCost::set_f_ref(const std::vector<Eigen::Vector3d>& f_ref) {
  f_ref_ = f_ref;
}


void ImpulseForceCost::set_f_weight(const std::vector<Eigen::Vector3d>& f_weight) {
  f_weight_ = f_weight;
}


double ImpulseForceCost::l(Robot& robot, const ContactStatus& contact_status, 
                           CostFunctionData& data, const double t, 
                           const ImpulseSplitSolution& s) const {
  double l = 0;
  for (int i=0; i<s.f.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      l+= (f_weight_[i].array()*(s.f[i]-f_ref_[i]).array()*(s.f[i]-f_ref_[i]).array()).sum();
    }
  }
  return 0.5 * l;
}


void ImpulseForceCost::lf(Robot& robot, const ContactStatus& contact_status, 
                          CostFunctionData& data, const double t, 
                          const std::vector<Eigen::Vector3d>& f, 
                          Eigen::VectorXd& lf) const {
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      lf.template segment<3>(dimf_stack).array() 
          += (f_weight_[i].array()*(f[i]-f_ref_[i]).array());
      dimf_stack += 3;
    }
  }
}


void ImpulseForceCost::lff(Robot& robot, const ContactStatus& contact_status, 
                           CostFunctionData& data, const double t, 
                           const std::vector<Eigen::Vector3d>& f, 
                           Eigen::MatrixXd& Qff) const {
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Qff.diagonal().template segment<3>(dimf_stack).noalias() += f_weight_[i];
      dimf_stack += 3;
    }
  }
}

} // namespace idocp