#include "robotoc/cost/trot_swing_foot_ref.hpp"

#include <cmath>


namespace robotoc {

TrotSwingFootRef::TrotSwingFootRef(const int contact_index,
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


TrotSwingFootRef::TrotSwingFootRef(const Robot& robot, 
                                   const std::string& contact_frame_name,
                                   const int x_ref_foot_contact_index, 
                                   const int y_ref_foot_contact_index, 
                                   const double step_length, 
                                   const double step_height) 
  : TrotSwingFootRef(robot.createContactStatus().findContactIndex(contact_frame_name),
                     x_ref_foot_contact_index, y_ref_foot_contact_index,
                     step_length, step_height) {
}


TrotSwingFootRef::~TrotSwingFootRef() {
}


void TrotSwingFootRef::updateRef(const ContactStatus& contact_status, 
                                      Eigen::VectorXd& x3d_ref) const {
  constexpr double eps = std::numeric_limits<double>::epsilon();
  const double xdiff = contact_status.contactPosition(contact_index_).coeff(0) 
                        - contact_status.contactPosition(x_ref_foot_contact_index_).coeff(0);
  if (std::abs(xdiff) < eps) {
    // This means that the foot step is the first step and the step length is half.
    x3d_ref.coeffRef(0) = contact_status.contactPosition(x_ref_foot_contact_index_).coeff(0);
    x3d_ref.coeffRef(0) += 0.25 * step_length_;
  }
  else {
    x3d_ref.coeffRef(0) = contact_status.contactPosition(x_ref_foot_contact_index_).coeff(0);
  }
  x3d_ref.coeffRef(1) = contact_status.contactPosition(y_ref_foot_contact_index_).coeff(1);
  x3d_ref.coeffRef(2) = step_height_;
}

} // namespace robotoc