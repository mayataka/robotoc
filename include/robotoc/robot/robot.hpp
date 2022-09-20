#ifndef ROBOTOC_ROBOT_HPP_
#define ROBOTOC_ROBOT_HPP_

#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/spatial/force.hpp"

#include "robotoc/robot/robot_model_info.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/robot/point_contact.hpp"
#include "robotoc/robot/surface_contact.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/robot/robot_properties.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

///
/// @class Robot
/// @brief Dynamics and kinematics model of robots. Wraps pinocchio::Model and 
/// pinocchio::Data. Includes contacts.
///
class Robot {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Constructs a robot model according to the input model info.
  /// @param[in] robot_model_info Info of a robot model.
  ///
  Robot(const RobotModelInfo& robot_model_info);

  ///
  /// @brief Default constructor. 
  ///
  Robot();

  ///
  /// @brief Default destructor. 
  ///
  ~Robot() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  Robot(const Robot&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  Robot& operator=(const Robot&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  Robot(Robot&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  Robot& operator=(Robot&&) noexcept = default;

  ///
  /// @brief Integrates the generalized velocity, that is, performs
  /// q <- q + integration_length * v like computation on the configuration 
  /// manifolds.
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] integration_length The length of the integration.
  /// @param[in, out] q Configuration. Size must be Robot::dimq().
  ///
  template <typename TangentVectorType, typename ConfigVectorType>
  void integrateConfiguration(
      const Eigen::MatrixBase<TangentVectorType>& v, 
      const double integration_length, 
      const Eigen::MatrixBase<ConfigVectorType>& q) const;

  ///
  /// @brief Integrates the generalized velocity, that is, performs
  /// q_integrated = q + integration_length * v like computation on the 
  /// configuration manifolds.
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] integration_length The length of the integration.
  /// @param[out] q_integrated Configuration. Size must be Robot::dimq().
  ///
  template <typename ConfigVectorType1, typename TangentVectorType,  
            typename ConfigVectorType2>
  void integrateConfiguration(
      const Eigen::MatrixBase<ConfigVectorType1>& q, 
      const Eigen::MatrixBase<TangentVectorType>& v, 
      const double integration_length, 
      const Eigen::MatrixBase<ConfigVectorType2>& q_integrated) const;

  ///
  /// @brief Transport the Jacobian w.r.t. the configuration evaluated at 
  /// the integrated configuration q + v to q. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] Jin The Jacobian w.r.t. the configuration evaluated at 
  /// q + v. Cols must be Robot::dimv().
  /// @param[out] Jout The resultant Jacobian w.r.t. the configuration evaluated 
  /// at q. The number of columns must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType, 
            typename MatrixType1, typename MatrixType2>
  void dIntegrateTransport_dq(
      const Eigen::MatrixBase<ConfigVectorType>& q,
      const Eigen::MatrixBase<TangentVectorType>& v,
      const Eigen::MatrixBase<MatrixType1>& Jin,
      const Eigen::MatrixBase<MatrixType2>& Jout) const;

  ///
  /// @brief Transport the Jacobian w.r.t. the generalized velocity evaluated at 
  /// the integrated configuration q + v to q. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] Jin The Jacobian w.r.t. the generalized velocity evaluated at 
  /// q + v. Cols must be Robot::dimv().
  /// @param[out] Jout The resultant Jacobian w.r.t. the generalized velocity 
  /// evaluated at q. The number of columns must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType, 
            typename MatrixType1, typename MatrixType2>
  void dIntegrateTransport_dv(
      const Eigen::MatrixBase<ConfigVectorType>& q,
      const Eigen::MatrixBase<TangentVectorType>& v,
      const Eigen::MatrixBase<MatrixType1>& Jin,
      const Eigen::MatrixBase<MatrixType2>& Jout) const;

  ///
  /// @brief Computes qf - q0 at the tangent space. 
  /// The result means that the unit velocity from initial configuration
  /// q0 to terminal configuration qf.
  /// @param[in] qf Terminal configuration. Size must be Robot::dimq().
  /// @param[in] q0 Initial configuration. Size must be Robot::dimq().
  /// @param[out] qdiff Difference of the configurations, qf - q0. Size must be 
  /// Robot::dimv().
  ///
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename TangentVectorType>
  void subtractConfiguration(
      const Eigen::MatrixBase<ConfigVectorType1>& qf, 
      const Eigen::MatrixBase<ConfigVectorType2>& q0,
      const Eigen::MatrixBase<TangentVectorType>& qdiff) const;

  ///
  /// @brief Computes the partial derivative of the function of subtraction
  /// qf - q0 w.r.t. the terminal configuration qf. 
  /// @param[in] qf Terminal configuration. Size must be Robot::dimq().
  /// @param[in] q0 Initial configuration. Size must be Robot::dimq().
  /// @param[out] dqdiff_dqf The partial derivative of the subtraction.
  /// Size must be Robot::dimv() x Robot::dimv().
  ///
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename MatrixType>
  void dSubtractConfiguration_dqf(
      const Eigen::MatrixBase<ConfigVectorType1>& qf,
      const Eigen::MatrixBase<ConfigVectorType2>& q0,
      const Eigen::MatrixBase<MatrixType>& dqdiff_dqf) const;

  ///
  /// @brief Computes the partial derivative of the function of subtraction
  /// qf - q0 w.r.t. the initial configuration qf. 
  /// @param[in] qf Terminal configuration. Size must be Robot::dimq().
  /// @param[in] q0 Initial configuration. Size must be Robot::dimq().
  /// @param[out] dqdiff_dq0 The partial derivative of the subtraction.
  /// Size must be Robot::dimv() x Robot::dimv().
  ///
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename MatrixType>
  void dSubtractConfiguration_dq0(
      const Eigen::MatrixBase<ConfigVectorType1>& qf,
      const Eigen::MatrixBase<ConfigVectorType2>& q0,
      const Eigen::MatrixBase<MatrixType>& dqdiff_dq0) const;

  ///
  /// @brief Computes an interpolation of two configurations by
  /// (1.0 - t) * s1 + t * s2
  /// @param[in] q1 Configuration. Size must be Robot::dimq().
  /// @param[in] q2 Configuration. Size must be Robot::dimq().
  /// @param[in] t Interpolation parameter. Must be between 0 to 1.
  /// @param[out] qout Configuration. Size must be Robot::dimq().
  ///
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename ConfigVectorType3>
  void interpolateConfiguration(
      const Eigen::MatrixBase<ConfigVectorType1>& q1, 
      const Eigen::MatrixBase<ConfigVectorType2>& q2, const double t,
      const Eigen::MatrixBase<ConfigVectorType3>& qout) const;

  ///
  /// @brief Transforms the Jacobian w.r.t. configuration to tangentspace.
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in, out] J Jacobian. Size must be Robot::dimq() x Robot::dimv().
  ///
  template <typename ConfigVectorType, typename MatrixType>
  void integrateCoeffWiseJacobian(const Eigen::MatrixBase<ConfigVectorType>& q, 
                                  const Eigen::MatrixBase<MatrixType>& J) const;

  ///
  /// @brief Updates the kinematics of the robot. The frame placements, frame 
  /// velocity, frame acceleration, and the relevant Jacobians are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] a Generalized acceleration. Size must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2>
  void updateKinematics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                        const Eigen::MatrixBase<TangentVectorType1>& v, 
                        const Eigen::MatrixBase<TangentVectorType2>& a);

  ///
  /// @brief Updates the kinematics of the robot. The frame placements, frame 
  /// velocity, and the relevant Jacobians are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType>
  void updateKinematics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                        const Eigen::MatrixBase<TangentVectorType>& v);

  ///
  /// @brief Updates the kinematics of the robot. The frame placements and
  /// and the relevant Jacobians are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  ///
  template <typename ConfigVectorType>
  void updateKinematics(const Eigen::MatrixBase<ConfigVectorType>& q);

  ///
  /// @brief Updates the frame kinematics of the robot. The frame placements, 
  /// velocities, and accelerations are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] a Generalized acceleration. Size must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2>
  void updateFrameKinematics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                             const Eigen::MatrixBase<TangentVectorType1>& v, 
                             const Eigen::MatrixBase<TangentVectorType2>& a);

  ///
  /// @brief Updates the frame kinematics of the robot. The frame placements 
  /// and velocities are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType>
  void updateFrameKinematics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                             const Eigen::MatrixBase<TangentVectorType>& v);

  ///
  /// @brief Updates the frame kinematics of the robot. The frame placements 
  /// are calculated. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  ///
  template <typename ConfigVectorType>
  void updateFrameKinematics(const Eigen::MatrixBase<ConfigVectorType>& q);

  ///
  /// @brief Returns the position of the frame. Before calling this function, 
  /// updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @return Const reference to the position of the frame.
  ///
  const Eigen::Vector3d& framePosition(const int frame_id) const;

  ///
  /// @brief Returns the position of the frame. Before calling this function, 
  /// updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @return Const reference to the position of the frame.
  ///
  const Eigen::Vector3d& framePosition(const std::string& frame_name) const;

  ///
  /// @brief Returns the rotation matrix of the frame. Before calling this  
  /// function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @return Const reference to the rotation matrix of the frame.
  ///
  const Eigen::Matrix3d& frameRotation(const int frame_id) const;

  ///
  /// @brief Returns the rotation matrix of the frame. Before calling this  
  /// function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @return Const reference to the rotation matrix of the frame.
  ///
  const Eigen::Matrix3d& frameRotation(const std::string& frame_name) const;

  ///
  /// @brief Returns the SE(3) of the frame. Before calling this function, 
  /// updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @return Const reference to the SE(3) of the frame.
  ///
  const SE3& framePlacement(const int frame_id) const;

  ///
  /// @brief Returns the SE(3) of the frame. Before calling this function, 
  /// updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @return Const reference to the SE(3) of the frame.
  ///
  const SE3& framePlacement(const std::string& frame_name) const;

  ///
  /// @brief Returns the position of the center of mass. Before calling this 
  /// function, updateKinematics() or updateFrameKinematics() must be called.
  /// 
  const Eigen::Vector3d& CoM() const;

  ///
  /// @brief Returns the linear (translational) velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Translational velocity of the frame.
  ///
  Eigen::Vector3d frameLinearVelocity(
      const int frame_id, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the linear (translational) velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Translational velocity of the frame.
  ///
  Eigen::Vector3d frameLinearVelocity(
      const std::string& frame_name, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the angular velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Angular velocity of the frame.
  ///
  Eigen::Vector3d frameAngularVelocity(
      const int frame_id, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the angular velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Angular velocity of the frame.
  ///
  Eigen::Vector3d frameAngularVelocity(
      const std::string& frame_name, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the spatial velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_id Index of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Spatial velocity of the frame.
  ///
  Vector6d frameSpatialVelocity(
      const int frame_id, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the spatial velocity of the frame. Before calling  
  /// this function, updateKinematics() or updateFrameKinematics() must be called.
  /// @param[in] frame_name Name of the frame.
  /// @param[in] reference_frame Reference frame. Default is 
  /// pinocchio::LOCAL_WORLD_ALIGNED (world velocity).
  /// @return Spatial velocity of the frame.
  ///
  Vector6d frameSpatialVelocity(
      const std::string& frame_name, 
      const pinocchio::ReferenceFrame reference_frame=pinocchio::LOCAL_WORLD_ALIGNED) const;

  ///
  /// @brief Returns the linear velocity of the center of mass expressed in the 
  /// world frame. Before calling this function, updateKinematics() or 
  /// updateFrameKinematics() must be called.
  /// 
  const Eigen::Vector3d& CoMVelocity() const;

  ///
  /// @brief Computes the Jacobian of the frame position expressed in the local 
  /// coordinate. Before calling this function, updateKinematics() must be 
  /// called.
  /// @param[in] frame_id Index of the frame.
  /// @param[out] J Jacobian. Size must be 6 x Robot::dimv().
  ///
  template <typename MatrixType>
  void getFrameJacobian(const int frame_id, 
                        const Eigen::MatrixBase<MatrixType>& J);

  ///
  /// @brief Gets the Jacobian of the position of the center of mass. Before 
  /// calling this function, updateKinematics() must be called.
  /// @param[out] J Jacobian. Size must be 3 x Robot::dimv().
  ///
  template <typename MatrixType>
  void getCoMJacobian(const Eigen::MatrixBase<MatrixType>& J) const;

  ///
  /// @brief Transforms 3D quantity from local frame to the world frame. 
  /// @param[in] frame_id Index of the frame whose local coordinate is of 
  /// interest.
  /// @param[in] vec_local 3-dimensional vector-valued quantity expressed at the 
  /// local frame.
  /// @param[out] vec_world 3-dimensional vector-valued quantity expressed at 
  /// the world frame.
  ///
  template <typename Vector3dType>
  void transformFromLocalToWorld(
      const int frame_id, const Eigen::Vector3d& vec_local, 
      const Eigen::MatrixBase<Vector3dType>& vec_world) const;

  ///
  /// @brief Jacobian of the transformation mapping of a 3D quantity from local 
  /// frame to the world frame. 
  /// @param[in] frame_id Index of the frame whose local coordinate is of 
  /// interest.
  /// @param[in] vec_world 3-dimensional vector-valued quantity expressed at 
  /// the world frame.
  /// @param[out] J The Jacobian. Size must be 6 x Robot::dimv(). The result is
  /// stored at the top 3 rows. 
  ///
  template <typename Vector3dType, typename MatrixType>
  void getJacobianTransformFromLocalToWorld(
      const int frame_id, const Eigen::MatrixBase<Vector3dType>& vec_world,
      const Eigen::MatrixBase<MatrixType>& J);

  ///
  /// @brief Computes the residual of the contact constriants represented by 
  /// Baumgarte's stabilization method. Before calling this function, 
  /// updateKinematics() must be called.
  /// @param[in] contact_status Contact status.
  /// @param[out] baumgarte_residual Residuals in the contact constraints. Size
  /// must be ContactStatus::dimf().
  ///
  template <typename VectorType>
  void computeBaumgarteResidual(
      const ContactStatus& contact_status, 
      const Eigen::MatrixBase<VectorType>& baumgarte_residual);

  ///
  /// @brief Computes the partial derivatives of the contact constriants 
  /// represented by the Baumgarte's stabilization method. 
  /// Before calling this function, updateKinematics() must be called. 
  /// @param[in] contact_status Contact status.
  /// @param[out] baumgarte_partial_dq The partial derivative  with respect to 
  /// the configuaration. Size must be ContactStatus::dimf() x Robot::dimv().
  /// @param[out] baumgarte_partial_dv The partial derivative  with respect to 
  /// the velocity. Size must be ContactStatus::dimf() x Robot::dimv().
  /// @param[out] baumgarte_partial_da The partial derivative  with respect to 
  /// the acceleration. Size must be ContactStatus::dimf() x be Robot::dimv().
  ///
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const ContactStatus& contact_status, 
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

  ///
  /// @brief Computes the residual of the impulse velocity constraint. Before 
  /// calling this function, updateKinematics() must be called.
  /// @param[in] impulse_status Impulse status.
  /// @param[out] velocity_residual Residuals in the contact velocity 
  /// constraint. Size must be ImpulseStatus::dimf().
  ///
  template <typename VectorType>
  void computeImpulseVelocityResidual(
      const ImpulseStatus& impulse_status, 
      const Eigen::MatrixBase<VectorType>& velocity_residual) const;

  ///
  /// @brief Computes the partial derivatives of the impulse velocity constraint.
  /// Before calling this function, updateKinematics() must be called. 
  /// @param[in] impulse_status Impulse status.
  /// @param[out] velocity_partial_dq The partial derivative with respect to the 
  /// configuaration. Size must be ImpulseStatus::dimf() x Robot::dimv(). 
  /// @param[out] velocity_partial_dv The partial derivative with respect to the 
  /// velocity. Size must be ImpulseStatus::dimf() x Robot::dimv().
  ///
  template <typename MatrixType1, typename MatrixType2>
  void computeImpulseVelocityDerivatives(
      const ImpulseStatus& impulse_status, 
      const Eigen::MatrixBase<MatrixType1>& velocity_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& velocity_partial_dv);

  ///
  /// @brief Computes the residual of the contact position constraint at the 
  /// impulse. Before calling this function, updateKinematics() must be called.
  /// @param[in] impulse_status Impulse status.
  /// @param[out] position_residual Residuals in the contact position constraint.
  /// Size must be ImpulseStatus::dimf().
  ///
  template <typename VectorType>
  void computeContactPositionResidual(
      const ImpulseStatus& impulse_status, 
      const Eigen::MatrixBase<VectorType>& position_residual);

  ///
  /// @brief Computes the partial derivative of the contact position at the 
  /// impulse. Before calling this  function, updateKinematics() must be called.
  /// @param[in] impulse_status Impulse status.
  /// @param[out] position_partial_dq The result of the partial derivative  
  /// with respect to the configuaration. Rows must be at least 3. Cols must 
  /// be Robot::dimv().
  ///
  template <typename MatrixType>
  void computeContactPositionDerivative(
      const ImpulseStatus& impulse_status, 
      const Eigen::MatrixBase<MatrixType>& position_partial_dq);

  ///
  /// @brief Set contact forces in this robot model for each active contacts.  
  /// @param[in] contact_status Contact status.
  /// @param[in] f The stack of the contact wrenches represented in the local 
  /// coordinate of the contact frame. Size must be Robot::maxNumContacts(). 
  /// 
  void setContactForces(const ContactStatus& contact_status, 
                        const std::vector<Vector6d>& f);

  ///
  /// @brief Set impulse forces in this robot model for each active impulses. 
  /// @param[in] impulse_status Impulse status.
  /// @param[in] f The stack of the impulse wrenches represented in the local 
  /// coordinate of the contact frame. Size must be Robot::maxNumContacts(). 
  /// 
  void setImpulseForces(const ImpulseStatus& impulse_status, 
                        const std::vector<Vector6d>& f);

  ///
  /// @brief Computes inverse dynamics, i.e., generalized torques corresponding 
  /// for given configuration, velocity, acceleration, and contact forces. 
  /// If the robot has contacts, update contact forces via setContactForces() 
  /// before calling this function.
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] a Generalized acceleration. Size must be Robot::dimv().
  /// @param[out] tau Generalized torques for fully actuated system. Size must 
  /// be Robot::dimv().
  ///
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3>
  void RNEA(const Eigen::MatrixBase<ConfigVectorType>& q, 
            const Eigen::MatrixBase<TangentVectorType1>& v, 
            const Eigen::MatrixBase<TangentVectorType2>& a, 
            const Eigen::MatrixBase<TangentVectorType3>& tau);

  ///
  /// @brief Computes the partial dervatives of the function of inverse dynamics 
  /// with respect to the configuration, velocity, and acceleration. If the 
  /// robot has contacts, update contact forces via setContactForces() before
  /// calling this function.
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Generalized velocity. Size must be Robot::dimv().
  /// @param[in] a Generalized acceleration. Size must be Robot::dimv().
  /// @param[out] dRNEA_partial_dq The partial derivative of inverse dynamics 
  /// with respect to the configuration. The size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] dRNEA_partial_dv The partial derivative of inverse dynamics 
  /// with respect to the velocity. The size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] dRNEA_partial_da The partial derivative of inverse dynamics 
  /// with respect to the acceleration. The size must be 
  /// Robot::dimv() x Robot::dimv().
  ///   
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename MatrixType1, 
            typename MatrixType2, typename MatrixType3>
  void RNEADerivatives(const Eigen::MatrixBase<ConfigVectorType>& q, 
                       const Eigen::MatrixBase<TangentVectorType1>& v, 
                       const Eigen::MatrixBase<TangentVectorType2>& a,
                       const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
                       const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_dv, 
                       const Eigen::MatrixBase<MatrixType3>& dRNEA_partial_da);

  ///
  /// @brief Computes the residual of the impulse dynamics for given 
  /// configuration and impulse change in the generalized velocity, and impulse 
  /// forces by using RNEA. Before call this function,   update contact forces 
  /// via setImpulseForces().
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] dv Impulse change in the velocity. Size must be Robot::dimv().
  /// @param[out] res Residual of impulse dynamics. Size must be Robot::dimv(). 
  ///
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2>
  void RNEAImpulse(const Eigen::MatrixBase<ConfigVectorType>& q, 
                   const Eigen::MatrixBase<TangentVectorType1>& dv,
                   const Eigen::MatrixBase<TangentVectorType2>& res);

  ///
  /// @brief Computes the partial dervatives of the function of impuse dynamics 
  /// with respect to configuration and impulse changes in the velocity. Before 
  /// calling this function, update contact forces via setImpulseForces().
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] dv Impulse change in the velocity. Size must be Robot::dimv().
  /// @param[out] dRNEA_partial_dq The partial derivative of impulse dynamics 
  /// with respect to the configuration. The size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] dRNEA_partial_ddv The partial derivative of impulse dynamics 
  /// with respect to the change in velocity. The size must be 
  /// Robot::dimv() x Robot::dimv().
  ///   
  template <typename ConfigVectorType, typename TangentVectorType, 
            typename MatrixType1, typename MatrixType2>
  void RNEAImpulseDerivatives(
      const Eigen::MatrixBase<ConfigVectorType>& q, 
      const Eigen::MatrixBase<TangentVectorType>& dv, 
      const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_ddv);

  ///
  /// @brief Computes the inverse of the joint inertia matrix M.
  /// @param[in] M Joint inertia matrix. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Minv Inverse of the joint inertia matrix M. Size must be 
  /// Robot::dimv() x Robot::dimv().
  ///   
  template <typename MatrixType1, typename MatrixType2>
  void computeMinv(const Eigen::MatrixBase<MatrixType1>& M, 
                   const Eigen::MatrixBase<MatrixType2>& Minv);

  ///
  /// @brief Computes the inverse of the contact dynamics matrix [[M J^T], [J O]].
  /// @param[in] M Joint inertia matrix. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[in] J Contact Jacobian. Size must be 
  /// ContactStatus::dimf() x Robot::dimv().
  /// @param[out] MJtJinv Inverse of the matrix [[M J^T], [J O]]. Size must be 
  /// (Robot::dimv() + ContactStatus::dimf()) x 
  /// (Robot::dimv() + ContactStatus::dimf()).
  ///   
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeMJtJinv(const Eigen::MatrixBase<MatrixType1>& M, 
                      const Eigen::MatrixBase<MatrixType2>& J,
                      const Eigen::MatrixBase<MatrixType3>& MJtJinv);

  ///
  /// @brief Generates feasible configuration randomly.
  /// @return The random and feasible configuration. Size is Robot::dimq().
  ///
  Eigen::VectorXd generateFeasibleConfiguration() const;

  ///
  /// @brief Normalizes a configuration vector.
  /// @param[in, out] q The normalized configuration. Size must be Robot::dimq().
  ///
  template <typename ConfigVectorType>
  void normalizeConfiguration(
      const Eigen::MatrixBase<ConfigVectorType>& q) const;

  ///
  /// @brief Gets the id of the specified frame.
  /// @param[in] frame_name Frame name of interest.
  /// @return id of the specified frame. Returns maximum number of the frames 
  /// if the frame whose name is frame_name does not exist.
  /// 
  int frameId(const std::string& frame_name) const;

  ///
  /// @brief Gets the name of the specified frame.
  /// @param[in] frame_id Frame id of interest.
  /// @return Name of the specified frame. 
  /// 
  std::string frameName(const int frame_id) const;

  ///
  /// @brief Returns the total mass of this robot model.
  /// @return The total mass of this robot model.
  ///
  double totalMass() const;

  ///
  /// @brief Returns the total weight of this robot model.
  /// @return The total weight of this robot model.
  ///
  double totalWeight() const;

  ///
  /// @brief Returns the dimensiton of the configuration.
  /// @return The dimensiton of the configuration.
  /// 
  int dimq() const;

  ///
  /// @brief Returns the dimensiton of the velocity, i.e, tangent space.
  /// @return The dimensiton of the velocity, i.e, the dimension of the tangent 
  /// space of the configuration space.
  /// 
  int dimv() const;

  ///
  /// @brief Returns the dimensiton of the velocity, i.e, tangent space.
  /// @return The dimensiton of the velocity, i.e, the dimension of the tangent 
  /// space of the configuration space.
  /// 
  int dimu() const;

  ///
  /// @brief Returns the maximum dimension of the contacts.
  /// @return The maximum dimension of the contacts.
  /// 
  int max_dimf() const;

  ///
  /// @brief Returns the dimensiton of the generalized torques corresponding to 
  /// the passive joints.
  /// @return The dimensiton of the generalized torques corresponding to the 
  //// passive joints.
  /// 
  int dim_passive() const;

  ///
  /// @brief Returns true if the robot has a floating base and false if not.
  /// @returns true if the robot has a floating base and false if not.
  /// 
  bool hasFloatingBase() const;

  ///
  /// @brief Returns the maximum number of the contacts.
  /// @return The maximum number of the contacts.
  /// 
  int maxNumContacts() const;

  ///
  /// @brief Returns the maximum number of the point contacts.
  /// @return The maximum number of the point contacts.
  /// 
  int maxNumPointContacts() const;

  ///
  /// @brief Returns the maximum number of the surface contacts.
  /// @return The maximum number of the surface contacts.
  /// 
  int maxNumSurfaceContacts() const;

  ///
  /// @brief Returns the type of the contact.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return Contact type. 
  ///
  ContactType contactType(const int contact_index) const;

  ///
  /// @brief Returns the types of the contacts.
  /// @return Contact types. 
  ///
  std::vector<ContactType> contactTypes() const;

  ///
  /// @brief Retruns the indices of the frames involving the point or surface 
  /// contacts. 
  /// @return Indices of the contact frames.
  /// 
  std::vector<int> contactFrames() const;

  ///
  /// @brief Retruns the names of the contact frames involving the point or 
  //// surface contacts. 
  /// @return Names of the contact frames.
  /// 
  std::vector<std::string> contactFrameNames() const;

  ///
  /// @brief Retruns the indices of the frames involving the contacts. 
  /// @return Indices of the contact frames.
  /// 
  std::vector<int> pointContactFrames() const;

  ///
  /// @brief Retruns the names of the frames involving the contacts. 
  /// @return Names of the contact frames.
  /// 
  std::vector<std::string> pointContactFrameNames() const;

  ///
  /// @brief Retruns the indices of the frames involving the surface contacts. 
  /// @return Indices of the surface contact frames.
  /// 
  std::vector<int> surfaceContactFrames() const;

  ///
  /// @brief Retruns the names of the frames involving the surface contacts. 
  /// @return Names of the surface contact frames.
  /// 
  std::vector<std::string> surfaceContactFrameNames() const;

  ///
  /// @brief Creates a ContactStatus for this robot model. 
  /// @return ContactStatus for this robot model.
  /// 
  ContactStatus createContactStatus() const;

  ///
  /// @brief Creates a ImpulseStatus for this robot model. 
  /// @return ImpulseStatus for this robot model.
  /// 
  ImpulseStatus createImpulseStatus() const;

  ///
  /// @brief Initializes the results of jointEffortLimit(), jointVelocityLimit(), 
  /// lowerJointPositionLimit(), and upperJointPositionLimit() by the URDF.
  /// 
  void initializeJointLimits();

  ///
  /// @brief Returns the effort limit of each joints.
  /// @return The effort limit of each joints.
  /// 
  Eigen::VectorXd jointEffortLimit() const;

  ///
  /// @brief Returns the joint velocity limit of each joints.
  /// @return The joint velocity limit of each joints.
  /// 
  Eigen::VectorXd jointVelocityLimit() const;

  ///
  /// @brief Returns the lower limit of the position of each joints.
  /// @return The lower limit of the position of each joints.
  /// 
  Eigen::VectorXd lowerJointPositionLimit() const;

  ///
  /// @brief Returns the upper limit of the position of each joints.
  /// @return The upper limit of the position of each joints.
  ///
  Eigen::VectorXd upperJointPositionLimit() const;

  ///
  /// @brief Sets the effort limit of each joints.
  /// @param[in] joint_effort_limit The effort limit of each joints.
  /// 
  void setJointEffortLimit(const Eigen::VectorXd& joint_effort_limit);

  ///
  /// @brief Sets the joint velocity limit of each joints.
  /// @param[in] joint_velocity_limit Sets the joint velocity limit of each 
  /// joints.
  ///
  void setJointVelocityLimit(const Eigen::VectorXd& joint_velocity_limit);

  ///
  /// @brief Sets the lower limit of the position of each joints.
  /// @param[in] lower_joint_position_limit The lower limit of the position of 
  /// each joints.
  /// 
  void setLowerJointPositionLimit(
      const Eigen::VectorXd& lower_joint_position_limit);

  ///
  /// @brief Sets the upper limit of the position of each joints.
  /// @param[in] upper_joint_position_limit The upper limit of the position of 
  /// each joints.
  ///
  void setUpperJointPositionLimit(
      const Eigen::VectorXd& upper_joint_position_limit);

  ///
  /// @brief Gets the robot model info. 
  /// @return const reference to a the robot model info.
  /// 
  const RobotModelInfo& robotModelInfo() const;

  ///
  /// @brief Gets a collectio of the properties of this robot model. 
  /// @return const reference to a collectio of the properties of this robot model.
  /// 
  const RobotProperties& robotProperties() const;

  ///
  ///
  /// @brief Sets a collection of the properties for this robot model. 
  /// @param[in] properties A collection of the properties for this robot model.
  ///
  void setRobotProperties(const RobotProperties& properties);

  ///
  /// @brief Displays the robot model onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const Robot& robot);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  // Robot model info
  RobotModelInfo info_;
  // Pinocchio models and datas, and runtime variables
  pinocchio::Model model_, impulse_model_;
  pinocchio::Data data_, impulse_data_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  Eigen::MatrixXd dimpulse_dv_; 
  // Contact models
  aligned_vector<PointContact> point_contacts_;
  aligned_vector<SurfaceContact> surface_contacts_;
  // Robot model variables
  int dimq_, dimv_, dimu_, dim_passive_, max_dimf_, max_num_contacts_;
  RobotProperties properties_;
  Eigen::VectorXd joint_effort_limit_, joint_velocity_limit_, 
                  lower_joint_position_limit_, upper_joint_position_limit_;
};

} // namespace robotoc

#include "robotoc/robot/robot.hxx"

#endif // ROBOTOC_ROBOT_HPP_ 