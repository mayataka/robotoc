#include "robotoc/cost/trotting_swing_foot_ref.hpp"


namespace robotoc {

TrottingSwingFootRef::TrottingSwingFootRef(const int x_ref_foot_contact_index, 
                                           const int y_ref_foot_contact_index, 
                                           const double step_length, 
                                           const double step_height) 
  : SwingFootRefBase(),
    x_ref_foot_contact_index_(x_ref_foot_contact_index), 
    y_ref_foot_contact_index_(y_ref_foot_contact_index),
    step_length_(step_length),
    step_height_(step_height) {
}


TrottingSwingFootRef::~TrottingSwingFootRef() {
}


void TrottingSwingFootRef::update_q_3d_ref(const ContactStatus& contact_status, 
                                           Eigen::VectorXd& q_3d_ref) const {
  q_3d_ref.coeffRef(0) = contact_status.contactPoint(x_ref_foot_contact_index_).coeff(0);
  q_3d_ref.coeffRef(0) += 0.5 * step_length_;
  q_3d_ref.coeffRef(1) = contact_status.contactPoint(y_ref_foot_contact_index_).coeff(1);
  q_3d_ref.coeffRef(2) = step_height_;
}

} // namespace robotoc