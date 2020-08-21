#include "idocp/cost/contact_cost.hpp"

#include <iostream>


namespace idocp {

ContactCost::ContactCost(const Robot& robot)
  : CostFunctionComponentBase(),
    max_dimf_(robot.max_dimf()),
    f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
    f_weight_(Eigen::VectorXd::Zero(robot.max_dimf())) {
}


ContactCost::ContactCost()
  : CostFunctionComponentBase(),
    max_dimf_(0),
    f_ref_(),
    f_weight_() {
}


ContactCost::~ContactCost() {
}


void ContactCost::set_f_ref(const Eigen::VectorXd& f_ref) {
  if (f_ref.size() == max_dimf_) {
    f_ref_ = f_ref;
  }
  else {
    std::cout << "invalid argment in set_f_ref(): size of f_ref must be " 
              << max_dimf_ << std::endl;
  }
}


void ContactCost::set_f_weight(const Eigen::VectorXd& f_weight) {
  if (f_weight.size() == max_dimf_) {
    f_weight_ = f_weight;
  }
  else {
    std::cout << "invalid argment in set_f_weight(): size of f_weight must be " 
              << max_dimf_ << std::endl;
  }
}


double ContactCost::l(const Robot& robot, CostFunctionData& data, 
                      const double t, const double dtau, 
                      const SplitSolution& s) const {
  double l = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        l += f_weight_.coeff(3*i+j) * (s.f.coeff(3*i+j)-f_ref_.coeff(3*i+j)) 
                                    * (s.f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
      }
    }
  }
  return 0.5 * dtau * l;
}


double ContactCost::phi(const Robot& robot, CostFunctionData& data, 
                        const double t, const SplitSolution& s) const {
  return 0;
}


void ContactCost::lf(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  int dimf = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        kkt_residual.lf().coeffRef(dimf+j) 
            += dtau * f_weight_.coeff(3*i+j) * (s.f.coeff(dimf+j)-f_ref_.coeff(3*i+j));
      }
      dimf += 3;
    }
  }
}


void ContactCost::lff(const Robot& robot, CostFunctionData& data, 
                      const double t, const double dtau, 
                      const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  int dimf = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        kkt_matrix.Qff().coeffRef(dimf+j, dimf+j) += dtau * f_weight_.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
}

} // namespace idocp