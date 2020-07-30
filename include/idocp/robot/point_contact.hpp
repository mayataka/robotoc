#ifndef IDOCP_POINT_CONTACT_HPP_
#define IDOCP_POINT_CONTACT_HPP_

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
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
  PointContact(const PointContact& other) = default;

  // Use default assign operator.
  PointContact& operator=(const PointContact& other) = default;

  // Resets the parameters of the Baumgarte's stabilization method.
  // Argments:
  //    baumgarte_weight_on_velocity : The weight parameter of the Baumgrate's 
  //      stabilization method on the velocity error. Must be positive.
  //    baumgarte_weight_on_position : The weight parameter of the Baumgrate's 
  //      stabilization method on the position error. Must be positive.
  void resetBaugrarteParameters(const double baumgarte_weight_on_velocity, 
                                const double baumgarte_weight_on_position);

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
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          Eigen::MatrixXd& Jacobian, 
                          const bool transpose=false);

  // Computes the 3xdimv contact Jacobian represented in the local coordinate 
  // of the contact frame. Before calling this function, you have to update the 
  // kinematics (with respect to the position, velocity, and acceleration) of 
  // the model in pinocchio::Data. The contact Jacobian is set in the block 
  // part of Jacobian.
  // Argments:
  //    model: Pinocchio model of the robot.
  //    data: Pinocchio data of the robot kinematics.
  //    block_begin_index: The initial index where the Jacobian is stored. 
  //      If transpose=true, the left side column index of the block matrix.
  //      If transpose=false, the top row index of the block matrix.
  //    Jacobian: Jacobian of the contact frame is stored in this matrix. 
  //      If transpose=true, rows must be dimv and cols must be at least 3. 
  //      If transpose=false, rows must be at leaste 3 and cols must be dimv. 
  //    transpose: flag for transposing the Jacobian or not. If true, the 
  //      Jacobian is transposed. If false, the Jacobian is not transposed, 
  //      i.e., the original Jacobian is returned.
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const int block_begin_index,
                          Eigen::MatrixXd& Jacobian,
                          const bool transpose=false);

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_residual: The vector result is stored in. Size must be 3.
  void computeBaumgarteResidual(const pinocchio::Model& model, 
                                const pinocchio::Data& data, 
                                Eigen::Vector3d& baumgarte_residual) const;

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data. The result is 
  // stored in baumgarte_residual with indices from result_begin until
  // result_begin+3.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    result_begin: The start index of the result.
  //    baumgarte_residual: The vector result is stored in. Size must be at 
  //      least 3.
  void computeBaumgarteResidual(const pinocchio::Model& model, 
                                const pinocchio::Data& data, 
                                const int result_begin,
                                Eigen::VectorXd& baumgarte_residual) const;

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
  void computeBaumgarteDerivatives(const pinocchio::Model& model, 
                                   pinocchio::Data& data, 
                                   Eigen::MatrixXd& baumgarte_partial_dq, 
                                   Eigen::MatrixXd& baumgarte_partial_dv, 
                                   Eigen::MatrixXd& baumgarte_partial_da);

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. Before calling this 
  // function, you have to update the kinematics of the model in 
  // pinocchio::Data. The derivatives are set in the block rows 
  // baumgarte_partial_dq, baumgarte_partial_dv, and baumgarte_partial_da, with 
  // the indices of rows start with block_rows_begin.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_partial_dq: Partial of contact constraints with respect to 
  //      the generalized configuration. The rows must be at least 3 and the 
  //      cols must be dimv.
  //    baumgarte_partial_dv: Partial of contact constraints with respect to 
  //      the generalized velocity. The rows must be at least 3 and the cols
  //      must be dimv.
  //    baumgarte_partial_da: Partial of contact constraints with respect to 
  //      the generalized acceleration. The rows must be at least 3 and the cols
  //      must be dimv.
  void computeBaumgarteDerivatives(const pinocchio::Model& model, 
                                   pinocchio::Data& data, 
                                   const int block_rows_begin,
                                   Eigen::MatrixXd& baumgarte_partial_dq, 
                                   Eigen::MatrixXd& baumgarte_partial_dv, 
                                   Eigen::MatrixXd& baumgarte_partial_da);
 
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

  // Returns dimv, the dimensiton of the generalized velocity.
  int dimv() const;

  // Returns baumgarte_weight_on_velocity_, the weight parameter on the 
  // contact velocity in the Baumgarte's stabilization method.
  double baumgarte_weight_on_velocity() const;

  // Returns baumgarte_weight_on_position, the weight parameter on the 
  // contact position in the Baumgarte's stabilization method.
  double baumgarte_weight_on_position() const;
  
  // Returns the contact point.
  Eigen::Vector3d contact_point() const;

  // Returns SE(3) translation from the parent joint to the contact frame.
  pinocchio::SE3 jXf() const;

private:
  bool is_active_;
  int contact_frame_id_, parent_joint_id_, dimv_;
  double baumgarte_weight_on_velocity_, baumgarte_weight_on_position_;
  Eigen::Vector3d contact_point_;
  pinocchio::SE3 jXf_;
  pinocchio::SE3::ActionMatrixType fXj_;
  Eigen::MatrixXd J_frame_, joint_v_partial_dq_, joint_a_partial_dq_, 
                  joint_a_partial_dv_, joint_a_partial_da_;
};

} // namespace idocp


#endif // IDOCP_POINT_CONTACT_HPP_ 