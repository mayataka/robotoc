#ifndef ROBOTOC_UTILS_TRAJECTORY_VIEWER_HPP_
#define ROBOTOC_UTILS_TRAJECTORY_VIEWER_HPP_

#include <string>
#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "gepetto/viewer/corba/client.hh"

#include "robotoc/robot/robot.hpp"


namespace robotoc {

///
/// @class TrajectoryViewer
/// @brief Viewer of the solution trajectory of the optimal control problems
/// based on gepetto-viewer-corba and pinocchio-gepetto-viewer.
///
class TrajectoryViewer {
public:
  ///
  /// @brief Constructs the viewer. The path to the robot description packages 
  /// is constructed by assuming that the URDF file lies as 
  /// "*/description_pkd_dir/foo/path_to_urdf".
  /// @param[in] path_to_urdf Path to the URDF.
  /// @param[in] base_joint_type Type of the base joint. Choose from 
  /// BaseJointType::FixedBase or BaseJointType::FloatingBase. Default is 
  /// BaseJointType::FixedBase.
  ///
  TrajectoryViewer(const std::string& path_to_urdf,
                   const BaseJointType& base_joint_type=BaseJointType::FixedBase);

  ///
  /// @brief Constructs the viewer.
  /// @param[in] path_to_urdf Path to the URDF.
  /// @param[in] path_to_pkg Path to the robot description packages.
  /// @param[in] base_joint_type Type of the base joint. Choose from 
  /// BaseJointType::FixedBase or BaseJointType::FloatingBase. Default is 
  /// BaseJointType::FixedBase.
  ///
  TrajectoryViewer(const std::string& path_to_urdf, 
                   const std::string& path_to_pkg,
                   const BaseJointType& base_joint_type=BaseJointType::FixedBase);

  ///
  /// @brief Displays the solution trajectory.
  /// @param[in] q_traj Std vector of the solution trajectory of the
  /// configuration q. (e.g., returned value of OCPSolver::getSolution("q") or
  /// ParNMPCSolver::getSolution("q")).
  /// @param[in] dt Time steps.
  ///
  void display(const std::vector<Eigen::VectorXd>& q_traj, 
               const std::vector<double>& dt);

  ///
  /// @brief Displays the solution trajectory.
  /// @param[in] q_traj Std vector of the solution trajectory of the
  /// configuration q. (e.g., returned value of OCPSolver::getSolution("q") or
  /// ParNMPCSolver::getSolution("q")).
  /// @param[in] dt Time step.
  /// @note All time steps are assumed to be equivalent.
  ///
  void display(const std::vector<Eigen::VectorXd>& q_traj, 
               const double dt);

  ///
  /// @brief Displays the solution trajectory with contact forces and friction
  /// cones (optional).
  /// @param[in] robot Robot model.
  /// @param[in] q_traj Std vector of the solution trajectory of the
  /// configuration q. (returned value of OCPSolver::getSolution("q") or
  /// ParNMPCSolver::getSolution("q"))
  /// @param[in] f_traj Std vector of the solution trajectory of the stack of 
  /// the contact forces expressed in the world frame f. If the contact is 
  /// inactive, the corresponding contact force must be set as 3D zero-vector.
  /// (e.g., returned value of OCPSolver::getSolution("f", "WORLD") or
  /// ParNMPCSolver::getSolution("f", "WORLD")).
  /// @param[in] dt Time steps.
  /// @param[in] mu Friction coefficient. If set by a positive value, the 
  /// friction cones are also displayed. If set by zero, they are not displayed.
  /// Default is zero.
  ///
  void display(Robot& robot, const std::vector<Eigen::VectorXd>& q_traj, 
               const std::vector<Eigen::VectorXd>& f_traj, 
               const std::vector<double>& dt, const double mu=0);

  ///
  /// @brief Displays the solution trajectory with contact forces and friction
  /// cones (optional).
  /// @param[in] robot Robot model.
  /// @param[in] q_traj Std vector of the solution trajectory of the
  /// configuration q. (returned value of OCPSolver::getSolution("q") or
  /// ParNMPCSolver::getSolution("q"))
  /// @param[in] f_traj Std vector of the solution trajectory of the stack of 
  /// the contact forces expressed in the world frame f. If the contact is 
  /// inactive, the corresponding contact force must be set as 3D zero-vector.
  /// (e.g., returned value of OCPSolver::getSolution("f", "WORLD") or
  /// ParNMPCSolver::getSolution("f", "WORLD")).
  /// @param[in] dt Time step.
  /// @param[in] mu Friction coefficient. If set by a positive value, the 
  /// friction cones are also displayed. If set by zero, they are not displayed.
  /// Default is zero.
  ///
  void display(Robot& robot, const std::vector<Eigen::VectorXd>& q_traj, 
               const std::vector<Eigen::VectorXd>& f_traj, const double dt, 
               const double mu=0);

  ///
  /// @brief Sets the camera transform. 
  /// @param[in] position Position of the camera. 
  /// @param[in] quat Quaternion of rotation of the camera.
  ///
  void setCameraTransform(const Eigen::Vector3d& position,
                          const Eigen::Vector4d& quat);

  ///
  /// @brief Sets the camera transform to the default position. 
  ///
  void setCameraTransformDefault();

  ///
  /// @brief Prints the current camera transform. 
  ///
  void printCurrentCameraTransform() const;

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
  pinocchio::Model model_;
  pinocchio::GeometryModel vmodel_;
  double force_radius_, force_length_, force_scale_, friction_cone_scale_;
  Eigen::Vector3d x_axis_, camera_pos_;
  Eigen::Vector4d camera_quat_;
  gepetto::corbaserver::Color force_color_, cone_color_;

  void setForceProperties();
  void setFrictionConeProperties();
  void setCameraTransform();

};

} // namespace robotoc

#endif // ROBOTOC_UTILS_TRAJECTORY_VIEWER_HPP_ 