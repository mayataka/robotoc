#include "idocp/utils/trajectory_viewer.hpp"


namespace idocp {

TrajectoryViewer::TrajectoryViewer(
    const std::string& path_to_raisim_activation_key,
    const std::string& path_to_urdf_for_raisim) 
  : path_to_raisim_activation_key_(path_to_raisim_activation_key),
    path_to_urdf_for_raisim_(path_to_urdf_for_raisim) {
}

} // namespace idocp