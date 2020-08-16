#ifndef IDOCP_POINT_CONTACT_HPP_
#define IDOCP_POINT_CONTACT_HPP_

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/spatial/force.hpp"


namespace idocp {

class PointContact {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Allocate matrices and vectors.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized, e.g., by pinocchio::buildModel() or 
  //      pinocchio::buildModelFromXML().
  //    contact_frame_id: The index of the contact frame. 
  //    baumgarte_alpha: The weight parameter of the Baumgrate's stabilization
  //      method on the velocity error. Must be positive.
  //    baumgarte_beta: The weight parameter of the Baumgrate's stabilization
  //      method on the position error. Must be positive.
  PointContact(const pinocchio::Model& model, const int contact_frame_id, 
               const double baumgarte_weight_on_velocity, 
               const double baumgarte_weight_on_position);

  // Default constructor.
  PointContact();

  // Destructor. 
  ~PointContact();

  // Use default copy constructor.
  PointContact(const PointContact&) = default;

  // Use default assign operator.
  PointContact& operator=(const PointContact&) = default;

  // Use default mvoe constructor.
  PointContact(PointContact&&) noexcept = default;

  // Use default move assign operator.
  PointContact& operator=(PointContact&&) noexcept = default;

  // Resets the parameters of the Baumgarte's stabilization method.
  // Argments:
  //    baumgarte_weight_on_velocity : The weight parameter of the Baumgrate's 
  //      stabilization method on the velocity error. Must be positive.
  //    baumgarte_weight_on_position : The weight parameter of the Baumgrate's 
  //      stabilization method on the position error. Must be positive.
  void resetBaugrarteParameters(const double baumgarte_weight_on_velocity, 
                                const double baumgarte_weight_on_position);

  // Resets the contact point.
  // Argments:
  //    contact_point: The contact point.
  void resetContactPoint(const Eigen::Vector3d& contact_point);

  // Resets the contact point by current kinematics of the robot. The kinematics
  // is passed through pinocchio::Data. Before calling this function, you have 
  // to update the kinematics (only with respect to the position) in 
  // pinocchio::Data.
  // Argments:
  //    data: The data including kinematics of the robot.
  void resetContactPointByCurrentKinematics(const pinocchio::Data& data);

  // Converts the contact forces to the corresponding joint forces.
  // Argments:
  //    contact_force: The contact forces in the world frame.
  //    joint_force: The joint forces of the robot model. 
  void computeJointForceFromContactForce(
      const Eigen::Vector3d& contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const;

  // Computes the 3xdimv contact Jacobian represented in the local coordinate 
  // of the contact frame. Before calling this function, you have to update the 
  // kinematics (with respect to the position, velocity, and acceleration) of 
  // the model in pinocchio::Data.
  // Argments:
  //    model: Pinocchio model of the robot.
  //    data: Pinocchio data of the robot kinematics.
  //    Jacobian: Jacobian of the contact frame is stored in this variable. 
  //      If transpose=true, size must be dimvx3. 
  //      If transpose=false, size must be 3xdimv. 
  //    transpose: flag for transposing the Jacobian or not. If true, the 
  //      Jacobian is transposed. If false, the Jacobian is not transposed, 
  //      i.e., the original Jacobian is returned.
  template <typename MatrixType>
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const Eigen::MatrixBase<MatrixType>& Jacobian,
                          const bool transpose=false);

  // Computes the 3xdimv contact Jacobian represented in the local coordinate 
  // of the contact frame. Before calling this function, you have to update the 
  // kinematics (with respect to the position, velocity, and acceleration) of 
  // the model in pinocchio::Data.
  // Argments:
  //    model: Pinocchio model of the robot.
  //    data: Pinocchio data of the robot kinematics.
  //    Jacobian: Jacobian of the contact frame is stored in this variable. 
  //      If transpose=true, size must be dimvx3. 
  //      If transpose=false, size must be 3xdimv. 
  //    transpose: flag for transposing the Jacobian or not. If true, the 
  //      Jacobian is transposed. If false, the Jacobian is not transposed, 
  //      i.e., the original Jacobian is returned.
  template <typename MatrixType>
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const double coeff, 
                          const Eigen::MatrixBase<MatrixType>& Jacobian,
                          const bool transpose=false);

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_residual: The vector result is stored in. Size must be 3.
  template <typename VectorType>
  void computeBaumgarteResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_residual: The vector result is stored in. Size must be 3.
  template <typename VectorType>
  void computeBaumgarteResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const double coeff, 
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. Before calling this 
  // function, you have to update the kinematics of the model in 
  // pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_partial_dq: Partial of contact constraints with respect to 
  //      the generalized configuration. size must be 3xdimv.
  //    baumgarte_partial_dv: Partial of contact constraints with respect to 
  //      the generalized velocity. size must be 3xdimv.
  //    baumgarte_partial_da: Partial of contact constraints with respect to 
  //      the generalized acceleration. size must be 3xdimv.
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const pinocchio::Model& model, pinocchio::Data& data, 
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. Before calling this 
  // function, you have to update the kinematics of the model in 
  // pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_partial_dq: Partial of contact constraints with respect to 
  //      the generalized configuration. size must be 3xdimv.
  //    baumgarte_partial_dv: Partial of contact constraints with respect to 
  //      the generalized velocity. size must be 3xdimv.
  //    baumgarte_partial_da: Partial of contact constraints with respect to 
  //      the generalized acceleration. size must be 3xdimv.
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const pinocchio::Model& model, pinocchio::Data& data, const double coeff,
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

  // Activate the contact.
  void activate();
  
  // Deactivate the contact.
  void deactivate();

  // Check if the contact is active or not. If the contact is active, return 
  // true. If the contact is not active, return false.
  bool isActive() const;

  // Returns contact_frame_id, the index of the contact frame.
  int contact_frame_id() const;

  // Returns parent_joint_id, the index of the parent joint of the contact 
  // frame.
  int parent_joint_id() const;

  // Returns baumgarte_weight_on_velocity_, the weight parameter on the 
  // contact velocity in the Baumgarte's stabilization method.
  double baumgarte_weight_on_velocity() const;

  // Returns baumgarte_weight_on_position, the weight parameter on the 
  // contact position in the Baumgarte's stabilization method.
  double baumgarte_weight_on_position() const;
  
  // Returns the contact point.
  Eigen::Vector3d contact_point() const;


private:
  bool is_active_;
  int contact_frame_id_, parent_joint_id_, dimv_;
  double baumgarte_weight_on_velocity_, baumgarte_weight_on_position_;
  Eigen::Vector3d contact_point_;
  pinocchio::SE3 jXf_;
  pinocchio::Motion v_frame_;
  Eigen::Matrix3d v_linear_skew_, v_angular_skew_;
  Eigen::Matrix<double, 6, Eigen::Dynamic> J_frame_, frame_v_partial_dq_, 
                                           frame_a_partial_dq_, 
                                           frame_a_partial_dv_, 
                                           frame_a_partial_da_;
};

} // namespace idocp

#include "idocp/robot/point_contact.hxx"

#endif // IDOCP_POINT_CONTACT_HPP_