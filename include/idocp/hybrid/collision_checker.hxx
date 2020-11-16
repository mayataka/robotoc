#ifndef IDOCP_COLLISION_CHECKER_HXX_ 
#define IDOCP_COLLISION_CHECKER_HXX_

#include "idocp/hybrid/collision_checker.hpp"

namespace idocp {

inline CollisionChecker::CollisionChecker(const Robot& robot)
  : contact_points_(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    contact_frames_indices_(robot.contactFramesIndices()) {
}


inline CollisionChecker::CollisionChecker()
  : contact_points_(),
    contact_frames_indices_() {
}


inline CollisionChecker::~CollisionChecker() {
}


inline std::vector<bool> CollisionChecker::check(Robot& robot, 
                                                 const Eigen::VectorXd& q) {
  robot.updateFrameKinematics(q);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    contact_points_[i] = robot.framePosition(contact_frames_indices_[i]);
  }
  std::vector<bool> is_impulse;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (contact_points_[i].coeff(2) <= 0) is_impulse.push_back(true);
    else is_impulse.push_back(false);
  }
  return is_impulse;
}


inline const std::vector<Eigen::Vector3d>& 
CollisionChecker::contactFramePositions() const {
  return contact_points_;
}


inline void CollisionChecker::printContactFramePositions() const {
  for (int i=0; i<contact_points_.size(); ++i) {
    std::cout << "contact index " << i << ": " 
              << contact_points_[i].transpose() << i << std::endl;
  }
}
  
} // namespace idocp

#endif // IDOCP_COLLISION_CHECKER_HXX_ 