#include "idocp/cost/contact_force_cost.hpp"

#include <iostream>


namespace idocp {

ContactForceCost::ContactForceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    max_point_contacts_(robot.max_point_contacts()),
    max_dimf_(robot.max_dimf()),
    f_ref_(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    f_weight_(robot.max_point_contacts(), Eigen::Vector3d::Zero()) {
}


ContactForceCost::ContactForceCost()
  : CostFunctionComponentBase(),
    max_point_contacts_(0),
    max_dimf_(0),
    f_ref_(),
    f_weight_() {
}


ContactForceCost::~ContactForceCost() {
}


bool ContactForceCost::useKinematics() const {
  return false;
}


void ContactForceCost::set_f_ref(const std::vector<Eigen::Vector3d>& f_ref) {
  if (f_ref.size() == max_point_contacts_) {
    f_ref_ = f_ref;
  }
  else {
    std::cout << "invalid argment in set_f_ref(): size of f_ref must be " 
              << max_point_contacts_ << std::endl;
  }
}


void ContactForceCost::set_f_weight(
    const std::vector<Eigen::Vector3d>& f_weight) {
  if (f_weight.size() == max_point_contacts_) {
    f_weight_ = f_weight;
  }
  else {
    std::cout << "invalid argment in set_f_weight(): size of f_ref must be " 
              << max_point_contacts_ << std::endl;
  }
}


double ContactForceCost::l(Robot& robot, CostFunctionData& data, const double t, 
                           const double dtau, const SplitSolution& s) const {
  double l = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      l += (f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array()) 
                                 * (s.f[i].array()-f_ref_[i].array())).sum();
    }
  }
  return 0.5 * dtau * l;
}


double ContactForceCost::phi(Robot& robot, CostFunctionData& data, 
                             const double t, const SplitSolution& s) const {
  return 0;
}


void ContactForceCost::lf(Robot& robot, CostFunctionData& data, const double t, 
                          const double dtau, const SplitSolution& s, 
                          SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += dtau * f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array());
      dimf_stack += 3;
    }
  }
}


void ContactForceCost::lff(Robot& robot, CostFunctionData& data, const double t, 
                           const double dtau, const SplitSolution& s, 
                           SplitKKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
          += dtau * f_weight_[i];
      dimf_stack += 3;
    }
  }
}

} // namespace idocp