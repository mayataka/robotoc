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
    max_dimf_(robot.max_dimf()),
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

} // namespace idocp