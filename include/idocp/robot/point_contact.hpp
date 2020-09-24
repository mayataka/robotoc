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

///
/// @class PointContact
/// @brief Kinematics model of a point contact.
///
class PointContact {
public:
  ///
  /// @brief Construct point contact model.
  /// @param[in] model The pinocchio model. Before call this constructor, 
  /// pinocchio model must be initialized, e.g., by pinocchio::buildModel().
  /// @param[in] contact_frame_id The index of the contact frame. 
  /// @param[in] friction_coefficient Friction coefficient. Must be positive.
  /// Defalut is 0.8.
  /// @param[in] restitution_coefficient Coefficient of the restitution. Must 
  /// be nonnegative. Must not be more than 1. Default is 0 
  /// (complete inelastic).
  ///
  PointContact(const pinocchio::Model& model, const int contact_frame_id, 
               const double friction_coefficient=0.8, 
               const double restitution_coefficient=0);

  ///
  /// @brief Default constructor. 
  ///
  PointContact();

  ///
  /// @brief Destructor. 
  ///
  ~PointContact();

  ///
  /// @brief Default copy constructor. 
  ///
  PointContact(const PointContact&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  PointContact& operator=(const PointContact&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  PointContact(PointContact&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  PointContact& operator=(PointContact&&) noexcept = default;

  ///
  /// @brief Converts the contact forces to the corresponding joint forces.
  /// @param[in] contact_force The contact forces in the local frame.
  /// @param[out] joint_forces: The corresponding joint forces of the robot 
  /// model. 
  ///
  void computeJointForceFromContactForce(
      const Eigen::Vector3d& contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const;

  ///
  /// @brief Computes the contact Jacobian represented in the local coordinate 
  /// of the contact frame. Before calling this function, you have to update the 
  /// kinematics (with respect to the position, velocity, and acceleration) of 
  /// the model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[out] Jacobian Jacobian of the contact frame is stored in this 
  /// variable. If transpose=true, size must be Robot::dimv() x 3. 
  /// If transpose=false, size must be 3 x Robot::dimv() 
  /// @param[in] transpose A flag for transposing the Jacobian or not. If true, 
  /// the Jacobian is transposed. If false, the Jacobian is not transposed, 
  /// i.e., the original Jacobian is returned. Default is false.
  ///
  template <typename MatrixType>
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const Eigen::MatrixBase<MatrixType>& Jacobian,
                          const bool transpose=false);

  ///
  /// @brief Computes the contact Jacobian represented in the local coordinate 
  /// of the contact frame multiplied by the user-defined constanct coefficient. 
  /// Before calling this function, you have to update the kinematics (with 
  /// respect to the position, velocity, and acceleration) of the model in 
  /// pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[in] coeff Coefficient multiplied to the resultant Jacobian.
  /// @param[out] Jacobian Jacobian of the contact frame is stored in this 
  /// variable. If transpose=true, size must be Robot::dimv() x 3. 
  /// If transpose=false, size must be 3 x Robot::dimv() 
  /// @param[in] transpose A flag for transposing the Jacobian or not. If true, 
  /// the Jacobian is transposed. If false, the Jacobian is not transposed, 
  /// i.e., the original Jacobian is returned. Default is false.
  ///
  template <typename MatrixType>
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const double coeff, 
                          const Eigen::MatrixBase<MatrixType>& Jacobian,
                          const bool transpose=false);

  ///
  /// @brief Computes the residual of the contact constraints considered by the 
  /// Baumgarte's stabilization method. Before calling this function, you have 
  /// to update the kinematics of the model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[in] time_step Time step of the discretization. Must be positive.
  /// @param[out] baumgarte_residual Residual of the Bamgarte's constraint.
  /// 
  template <typename VectorType>
  void computeBaumgarteResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const double time_step,
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  ///
  /// @brief Computes the residual of the contact constraints considered by the 
  /// Baumgarte's stabilization method multiplied by user-defined coefficient. 
  /// Before calling this function, you have to update the kinematics of the 
  /// model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[in] coeff Coefficient multiplied to the resultant Jacobian.
  /// @param[in] time_step Time step of the discretization. Must be positive.
  /// @param[out] baumgarte_residual Residual of the Bamgarte's constraint.
  /// 
  template <typename VectorType>
  void computeBaumgarteResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const double coeff, const double time_step,
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  ///
  /// @brief Computes the partial derivatives of the contact constraints
  /// considered by the Baumgarte's stabilization method. Before calling this 
  /// function, you have to update the kinematics of the model in 
  /// pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[in] time_step Time step of the discretization. Must be positive.
  /// @param[out] baumgarte_partial_dq The result of the partial derivative  
  /// with respect to the configuaration. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_dv The result of the partial derivative  
  /// with respect to the velocity. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_da The result of the partial derivative  
  /// with respect to the acceleration. Size must be 3 x Robot::dimv().
  /// 
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const pinocchio::Model& model, pinocchio::Data& data, 
      const double time_step,
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

  ///
  /// @brief Computes the partial derivatives of the contact constraints
  /// considered by the Baumgarte's stabilization method multiplied by 
  /// user-defined coefficient. Before calling this function, you have to 
  /// update the kinematics of the model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[in] coeff Coefficient multiplied to the resultant Jacobian.
  /// @param[in] time_step Time step of the discretization. Must be positive.
  /// @param[out] baumgarte_partial_dq The result of the partial derivative  
  /// with respect to the configuaration. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_dv The result of the partial derivative  
  /// with respect to the velocity. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_da The result of the partial derivative  
  /// with respect to the acceleration. Size must be 3 x Robot::dimv().
  /// 
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const pinocchio::Model& model, pinocchio::Data& data, const double coeff,
      const double time_step,
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

  ///
  /// @brief Computes the residual of the impulse velocity constraints 
  /// considered by the Baumgarte's stabilization method. Before calling this 
  /// function, you have to update the kinematics of the model in 
  /// pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[out] baumgarte_residual Residual of the Bamgarte's velocity 
  /// constraint.
  /// 
  template <typename VectorType>
  void computeBaumgarteImpulseResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  ///
  /// @brief Computes the partial derivatives of the contact velocity 
  /// constraints considered by the Baumgarte's stabilization method. 
  /// Before calling this function, you have to update the kinematics of the 
  /// model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[out] baumgarte_partial_dq The result of the partial derivative  
  /// with respect to the configuaration. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_dq The result of the partial derivative  
  /// with respect to the velocity. Size must be 3 x Robot::dimv().
  /// @param[out] baumgarte_partial_ddv The result of the partial derivative  
  /// with respect to the change in the velocity. Size must be 3 x Robot::dimv().
  /// 
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteImpulseDerivatives(
      const pinocchio::Model& model, pinocchio::Data& data,
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_ddv);

  ///
  /// @brief Computes the residual of the contact position constraints.
  /// Before calling this function, you have to update the kinematics of the  
  /// model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[out] contact_residual Residual of the contact constraint.
  /// 
  template <typename VectorType>
  void computeContactResidual(
      const pinocchio::Model& model, const pinocchio::Data& data, 
      const Eigen::MatrixBase<VectorType>& contact_residual) const;

  ///
  /// @brief Computes the partial derivative of the contact position  
  /// constraint with respect to the configuration. 
  /// Before calling this function, you have to update the kinematics of the 
  /// model in pinocchio::Data.
  /// @param[in] model Pinocchio model of the robot.
  /// @param[in] data Pinocchio data of the robot kinematics.
  /// @param[out] contact_partial_dq The result of the partial derivative  
  /// with respect to the configuaration. Size must be 3 x Robot::dimv().
  /// 
  template <typename MatrixType>
  void computeContactDerivative(
      const pinocchio::Model& model, pinocchio::Data& data,
      const Eigen::MatrixBase<MatrixType>& contact_partial_dq);

  ///
  /// @brief Activate the contact.
  ///
  void activate();
  
  ///
  /// @brief Deactivate the contact.
  ///
  void deactivate();

  ///
  /// @brief Check if the contact is active or not. If the contact is active, 
  /// return true. If the contact is not active, return false.
  /// @return true if the contact is active. false if the contact is not active.
  /// 
  bool isActive() const;

  ///
  /// @brief Sets the contact points.
  /// @param[in] contact_point The contact points, i.e., (x,y,z) position of 
  /// the end-effector frame having the contact.
  /// 
  void setContactPoint(const Eigen::Vector3d& contact_point);

  ///
  /// @brief Sets the contact points by current kinematics of the robot. 
  /// The kinematics is passed through pinocchio::Data. Before calling this 
  /// function, kinematics (only with respect to the position) in 
  /// pinocchio::Data must be updated.
  /// @param[in] data The data including kinematics of the robot.
  ///
  void setContactPointByCurrentKinematics(const pinocchio::Data& data);

  ///
  /// @brief Sets the friction coefficient.
  /// @param[in] friction_coefficient Friction coefficient. Must be positive.
  ///
  void setFrictionCoefficient(const double friction_coefficient);

  ///
  /// @brief Sets the coefficient of the resititution.
  /// @param[in] restitution_coefficient Coefficient of the restitution. Must 
  /// be nonnegative. Must not be more than 1. Default is 0 
  /// (complete inelastic).
  ///
  void setRestitutionCoefficient(const double restitution_coefficient);

  ///
  /// @brief Returns the contact point.
  /// @return Const reference to the contact point.
  ///
  const Eigen::Vector3d& contactPoint() const;

  ///
  /// @brief Returns friction coefficient.
  /// @return Friction coefficient.
  ///
  double frictionCoefficient() const;

  ///
  /// @brief Returns coefficient of restitution of Newton's impact law.
  /// @return Coefficient of restitution.
  ///
  double restitutionCoefficient() const;

  ///
  /// @brief Returns contact frame id, the index of the contact frame.
  /// @return Contact frame id.
  /// 
  int contact_frame_id() const;

  ///
  /// @brief Returns parent joint id, the index of the parent joint of the 
  /// contact frame.
  /// @return Parent joint id.
  /// 
  int parent_joint_id() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool is_active_;
  int contact_frame_id_, parent_joint_id_, dimv_;
  double friction_coefficient_, restitution_coefficient_;
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