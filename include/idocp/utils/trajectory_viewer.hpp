#ifndef IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_
#define IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"
#include "pinocchio/gepetto/viewer.hpp"
#include "gepetto/viewer/corba/client.hh"

#include "idocp/robot/robot.hpp"


namespace idocp {

class TrajectoryViewer {
public:
  TrajectoryViewer(const std::string& description_pkg_serach_path, 
                   const std::string& path_to_urdf);

  TrajectoryViewer(const std::string& path_to_urdf);

  void display(const std::vector<Eigen::VectorXd>& q_traj, 
               const double sampling_period_in_sec);

  void display(Robot& robot, const double mu, 
               const std::vector<Eigen::VectorXd>& q_traj, 
               const std::vector<Eigen::VectorXd>& f_traj,
               const double sampling_period_in_sec);

  static Eigen::Matrix3d rotationMatrix(const Eigen::Vector3d& a, 
                                        const Eigen::Vector3d& b) {
    const Eigen::Vector3d a_nrmrzd = a.normalized();
    const Eigen::Vector3d b_nrmrzd = b.normalized();
    const Eigen::Vector3d a_cross_b = a_nrmrzd.cross(b_nrmrzd);
    const double s = a_cross_b.norm();
    if (s == 0.) {
      return Eigen::Matrix3d::Identity();
    }
    else {
      const double c = a_nrmrzd.dot(b_nrmrzd);
      Eigen::Matrix3d ab_skew(Eigen::Matrix3d::Zero());
      pinocchio::skew(a_cross_b, ab_skew);
      return Eigen::Matrix3d::Identity() + ab_skew + ab_skew * ab_skew * ((1-c) / (s*s));
    }
  }

private:
  std::string description_pkg_serach_path_, path_to_urdf_;
  double force_radius_, force_length_, friction_cone_scale_;
  Eigen::Vector3d x_axis_;
  gepetto::corbaserver::Color cone_color_;
  
};

} // namespace idocp

#endif // IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_ 