#ifndef IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_
#define IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
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

  Eigen::Matrix3d rotationMatrix(const Eigen::Vector3d& vec) const {
    const Eigen::Vector3d vec_nrmrzd = vec.normalized();
    const Eigen::Vector3d axis_cross_vec = x_axis_.cross(vec_nrmrzd);
    const double cross_norm = axis_cross_vec.norm();
    if (cross_norm == 0.) {
      return Eigen::Matrix3d::Identity();
    }
    else {
      const double inner = x_axis_.dot(vec_nrmrzd);
      Eigen::Matrix3d skew_mat(Eigen::Matrix3d::Zero());
      pinocchio::skew(axis_cross_vec, skew_mat);
      return Eigen::Matrix3d::Identity() 
              + skew_mat 
              + skew_mat * skew_mat * ((1-inner) / (cross_norm*cross_norm));
    }
  }

private:
  std::string description_pkg_serach_path_, path_to_urdf_;
  pinocchio::Model model_;
  pinocchio::GeometryModel vmodel_;
  double force_radius_, force_length_, friction_cone_scale_;
  Eigen::Vector3d x_axis_;
  gepetto::corbaserver::Color cone_color_;

  void setForceProperties();
  void setFrictionConeProperties();
};

} // namespace idocp

#endif // IDOCP_UTILS_TRAJECTORY_VIEWER_HPP_ 