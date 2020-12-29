#ifndef IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_
#define IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_

#include "Eigen/Core"

#include "pinocchio/gepetto/viewer.hpp"

namespace idocp {

class TrajectoryViewer {
public:
  TrajectoryViewer(const std::string& path_to_urdf);

  void display(const std::vector<Eigen::VectorXd>& q_traj, 
               const double sampling_period_in_sec, const bool recording=false,
               const std::string& path_to_save_file="");

private:
  
};

} // namespace idocp

#endif // IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_ 