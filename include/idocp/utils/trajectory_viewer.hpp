#ifndef IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_
#define IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_

#include "Eigen/Core"


namespace idocp {

class TrajectoryViewer {
public:
  TrajectoryViewer(const std::string& path_to_raisim_activation_key,
                   const std::string& path_to_urdf_for_raisim);

  void display(const std::string& path_to_data, 
               const double sampling_period_in_sec, const bool recording=false,
               const std::string& path_to_save_file="");

  static void RaiSimOgreCallback();

private:
  std::string path_to_raisim_activation_key_, path_to_urdf_for_raisim_;

};

} // namespace idocp

#endif // IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_ 