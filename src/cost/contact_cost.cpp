#include "cost/contact_cost.hpp"

#include <assert.h>


namespace idocp {

ContactCost::ContactCost(const Robot& robot, const Eigen::VectorXd& f_weight)
  : contact_status_(robot.max_point_contacts(), false),
    max_point_contacts_(robot.max_point_contacts()),
    max_dimf_(robot.max_dimf()),
    f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
    f_weight_(f_weight) {
  assert(f_weight.size()==robot.max_dimf());
}


ContactCost::ContactCost(const Robot& robot, const Eigen::VectorXd& f_ref, 
                         const Eigen::VectorXd& f_weight)
  : contact_status_(robot.max_point_contacts(), false),
    max_point_contacts_(robot.max_point_contacts()),
    max_dimf_(robot.max_dimf()),
    f_ref_(f_ref),
    f_weight_(f_weight) {
  assert(f_ref.size()==robot.max_dimf());
  assert(f_weight.size()==robot.max_dimf());
}


void ContactCost::set_f_ref(const Eigen::VectorXd& f_ref) {
  assert(f_ref_.size()==max_dimf_);
  f_ref_ = f_ref;
}


void ContactCost::set_f_weight(const Eigen::VectorXd& f_weight) {
  assert(f_weight_.size()==max_dimf_);
  f_weight_ = f_weight;
}


void ContactCost::setContactStatus(const Robot& robot) {
  for (int i=0; i<max_point_contacts_; ++i) {
    contact_status_[i] = robot.is_contact_active(i);
  }
}


double ContactCost::l(const double dtau, const Eigen::VectorXd& f) const {
  assert(dtau > 0);
  assert(f.size()==max_dimf_);
  double l = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (contact_status_[i]) {
      for (int j=0; j<3; ++j) {
        l += f_weight_.coeff(3*i+j) * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j)) 
                                    * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
      }
    }
  }
  return 0.5 * dtau * l;
}


void ContactCost::lf(const double dtau, const Eigen::VectorXd& f, 
                     Eigen::VectorXd& lf) const {
  assert(dtau > 0);
  assert(f.size()==max_dimf_);
  assert(lf.size()==max_dimf_);
  int dimf = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (contact_status_[i]) {
      for (int j=0; j<3; ++j) {
        lf.coeffRef(dimf+j) = dtau * f_weight_.coeff(3*i+j) 
                                   * (f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
      }
      dimf += 3;
    }
  }
}


void ContactCost::lff(const double dtau, Eigen::MatrixXd& lff) const {
  assert(dtau > 0);
  assert(lff.rows() == max_dimf_);
  assert(lff.cols() == max_dimf_);
  int dimf = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (contact_status_[i]) {
      for (int j=0; j<3; ++j) {
        lff.coeffRef(dimf+j, dimf+j) = dtau * f_weight_.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
}


void ContactCost::augment_lff(const double dtau, Eigen::MatrixXd& lff) const {
  assert(dtau > 0);
  assert(lff.rows() == max_dimf_);
  assert(lff.cols() == max_dimf_);
  int dimf = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (contact_status_[i]) {
      for (int j=0; j<3; ++j) {
        lff.coeffRef(dimf+j, dimf+j) += dtau * f_weight_.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
}

} // namespace idocp