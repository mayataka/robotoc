#include "cost/contact_cost.hpp"

#include <assert.h>


namespace idocp {

ContactCost::ContactCost(const Robot& robot, const Eigen::VectorXd& f_weight)
  : f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
    f_weight_(f_weight) {
  assert(robot.max_dimf()==f_weight.size());
}


ContactCost::ContactCost(const Robot& robot, const Eigen::VectorXd& f_ref, 
                         const Eigen::VectorXd& f_weight)
  : f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
    f_weight_(f_weight) {
  assert(robot.max_dimf()==f_ref.size());
  assert(robot.max_dimf()==f_weight.size());
}


void ContactCost::set_f_ref(const Eigen::VectorXd& f_ref) {
  assert(f_ref_.size()==f_ref.size());
  f_ref_ = f_ref;
}


void ContactCost::set_f_weight(const Eigen::VectorXd& f_weight) {
  assert(f_weight_.size()==f_weight.size());
  f_weight_ = f_weight;
}


double ContactCost::l(const Robot& robot, const double dtau, 
                      const Eigen::VectorXd& f) {
  assert(dtau > 0);
  assert(f.size() == robot.max_dimf());
  double l = 0;
  for (int i=0; i<robot.dimf(); ++i) {
    l += f_weight_.coeff(i) * (f.coeff(i)-f_ref_.coeff(i)) 
                            * (f.coeff(i)-f_ref_.coeff(i));
  }
  return 0.5 * dtau * l;
}


void ContactCost::lf(const Robot& robot, const double dtau, 
                     const Eigen::VectorXd& f, Eigen::VectorXd& lf) {
  assert(dtau > 0);
  assert(f.size() == robot.max_dimf());
  assert(lf.size() == robot.max_dimf());
  for (int i=0; i<robot.dimf(); ++i) {
    lf.coeffRef(i) = dtau * f_weight_.coeff(i) * (f.coeff(i)-f_ref_.coeff(i));
  }
}


void ContactCost::lff(const Robot& robot, const double dtau, 
                      Eigen::MatrixXd& lff) {
  assert(dtau > 0);
  assert(lff.rows() == robot.max_dimf());
  assert(lff.cols() == robot.max_dimf());
  for (int i=0; i<robot.dimf(); ++i) {
    lff.coeffRef(i) += dtau * f_weight_.coeff(i);
  }
}

} // namespace idocp