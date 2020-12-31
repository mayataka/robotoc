#ifndef IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_
#define IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"
#include "pinocchio/gepetto/viewer.hpp"


namespace idocp {

class TrajectoryViewer {
public:
  TrajectoryViewer(const std::string& description_pkg_serach_path, 
                   const std::string& path_to_urdf);

  void display(const std::vector<Eigen::VectorXd>& q_traj, 
               const double sampling_period_in_sec);

private:
  std::string description_pkg_serach_path_, path_to_urdf_;
  
};

} // namespace idocp

#endif // IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_ 