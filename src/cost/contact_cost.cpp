#include "idocp/cost/contact_cost.hpp"

#include <iostream>


namespace idocp {

ContactCost::ContactCost(const Robot& robot)
  : CostFunctionComponentBase(),
    max_dimf_(robot.max_dimf()),
    f_ref_(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    f_weight_(robot.max_point_contacts(), Eigen::Vector3d::Zero()) {
}


ContactCost::ContactCost()
  : CostFunctionComponentBase(),
    max_dimf_(0),
    f_ref_(),
    f_weight_() {
}


ContactCost::~ContactCost() {
}


bool ContactCost::useKinematics() const {
  return false;
}


void ContactCost::set_f_ref(const std::vector<Eigen::Vector3d>& f_ref) {
  f_ref_ = f_ref;
}


void ContactCost::set_f_weight(const std::vector<Eigen::Vector3d>& f_weight) {
  f_weight_ = f_weight;
}


double ContactCost::l(Robot& robot, CostFunctionData& data, const double t, 
                      const double dtau, const SplitSolution& s) const {
  double l = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      l += (f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array()) 
                                 * (s.f[i].array()-f_ref_[i].array())).sum();
    }
  }
  return 0.5 * dtau * l;
}


double ContactCost::phi(Robot& robot, CostFunctionData& data, 
                        const double t, const SplitSolution& s) const {
  return 0;
}


void ContactCost::lf(Robot& robot, CostFunctionData& data, const double t, 
                     const double dtau, const SplitSolution& s, 
                     KKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      kkt_residual.lf().template segment<3>(dimf_stack).array()
          += dtau * f_weight_[i].array() * (s.f[i].array()-f_ref_[i].array());
      dimf_stack += 3;
    }
  }
}


void ContactCost::lff(Robot& robot, CostFunctionData& data, const double t, 
                      const double dtau, const SplitSolution& s, 
                      KKTMatrix& kkt_matrix) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
        kkt_matrix.Qff().diagonal().template segment<3>(dimf_stack).noalias() 
            += dtau * f_weight_[i];
      dimf_stack += 3;
    }
  }
}

} // namespace idocp