#include "robotoc/cost/trotting_swing_foot_ref.hpp"

#include <cmath>


namespace robotoc {

TrottingSwingFootRef::TrottingSwingFootRef(const int contact_index,
                                           const int x_ref_foot_contact_index, 
                                           const int y_ref_foot_contact_index, 
                                           const double step_length, 
                                           const double step_height) 
  : SwingFootRefBase(),
    contact_index_(contact_index),
    x_ref_foot_contact_index_(x_ref_foot_contact_index), 
    y_ref_foot_contact_index_(y_ref_foot_contact_index),
    step_length_(step_length),
    step_height_(step_height) {
}


TrottingSwingFootRef::~TrottingSwingFootRef() {
}


void TrottingSwingFootRef::update_q_3d_ref(const ContactStatus& contact_status, 
                                           Eigen::VectorXd& q_3d_ref) const {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  const double xdiff = contact_status.contactPoint(contact_index_).coeff(0) 
                        - contact_status.contactPoint(x_ref_foot_contact_index_).coeff(0);
  if (std::abs(xdiff) < eps) {
    // This means that the foot step is the first step and the step length is half.
    q_3d_ref.coeffRef(0) = contact_status.contactPoint(x_ref_foot_contact_index_).coeff(0);
    q_3d_ref.coeffRef(0) += 0.25 * step_length_;
  }
  else {
    q_3d_ref.coeffRef(0) = contact_status.contactPoint(x_ref_foot_contact_index_).coeff(0);
  }
  q_3d_ref.coeffRef(1) = contact_status.contactPoint(y_ref_foot_contact_index_).coeff(1);
  q_3d_ref.coeffRef(2) = step_height_;
}

} // namespace robotoc