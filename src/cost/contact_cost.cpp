#include "idocp/cost/contact_cost.hpp"

#include <assert.h>


namespace idocp {

ContactCost::ContactCost(const Robot& robot)
  : CostFunctionComponentBase(),
    f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
    f_weight_(Eigen::VectorXd::Zero(robot.max_dimf())) {
}


ContactCost::ContactCost()
  : f_ref_(),
    f_weight_() {
}


ContactCost::~ContactCost() {
}


void ContactCost::set_f_ref(const Eigen::VectorXd& f_ref) {
  assert(f_ref.size() == f_ref_.size());
  f_ref_ = f_ref;
}


void ContactCost::set_f_weight(const Eigen::VectorXd& f_weight) {
  assert(f_weight.size() == f_weight_.size());
  f_weight_ = f_weight;
}


double ContactCost::l(const Robot& robot, CostFunctionData& data, 
                      const double t, const double dtau, 
                      const Eigen::Ref<const Eigen::VectorXd> q, 
                      const Eigen::Ref<const Eigen::VectorXd> v, 
                      const Eigen::Ref<const Eigen::VectorXd> a, 
                      const Eigen::Ref<const Eigen::VectorXd> f, 
                      const Eigen::Ref<const Eigen::VectorXd> u) const {
  double l = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        l += f_weight_.coeff(3*i+j) * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j)) 
                                    * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
      }
    }
  }
  return 0.5 * dtau * l;
}


void ContactCost::lf(const Robot& robot, CostFunctionData& data, const double t, 
                     const double dtau, 
                     const Eigen::Ref<const Eigen::VectorXd> f, 
                     Eigen::Ref<Eigen::VectorXd> lf) const {
  int dimf = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lf.coeffRef(dimf+j) = dtau * f_weight_.coeff(3*i+j) 
                                   * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
      }
      dimf += 3;
    }
  }
}


void ContactCost::lff(const Robot& robot, CostFunctionData& data, 
                      const double t, const double dtau, 
                      const Eigen::Ref<const Eigen::VectorXd> f, 
                      Eigen::Ref<Eigen::MatrixXd> lff) const {
  assert(dtau > 0);
  assert(lff.rows() == robot.max_dimf());
  assert(lff.cols() == robot.max_dimf());
  int dimf = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lff.coeffRef(dimf+j, dimf+j) = dtau * f_weight_.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
}


void ContactCost::augment_lff(const Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const Eigen::Ref<const Eigen::VectorXd> f, 
                              Eigen::Ref<Eigen::MatrixXd> lff) const {
  assert(dtau > 0);
  assert(lff.rows() == robot.max_dimf());
  assert(lff.cols() == robot.max_dimf());
  int dimf = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lff.coeffRef(dimf+j, dimf+j) += dtau * f_weight_.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
}

} // namespace idocp