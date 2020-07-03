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

  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  //    contact_frame_id: The frame index of the contact. 
  //    baumgarte_alpha: The weight parameter of the Baumgrate's stabilization
  //      method
  //    baumgarte_beta: The weight parameter of the Baumgrate's stabilization
  //      method
  PointContact(const pinocchio::Model& model, 
               const unsigned int contact_frame_id, 
               const double baumgarte_alpha, const double baumgarte_beta);

  // Destructor. 
  ~PointContact();

  // Use default copy constructor.
  PointContact(const PointContact& other) = default;

  // Use default assign operator.
  PointContact& operator=(const PointContact& other) = default;

  // Move constructor.
  PointContact(PointContact&& other) noexcept;

  // Move assignment operator. 
  PointContact& operator=(PointContact&& other) noexcept;

  // Resets the parameters for the Baumgarte's stabilization method.
  // Argments:
  //    alpha: The weight parameter of the Baumgrate's stabilization method
  //    beta: The weight parameter of the Baumgrate's stabilization method
  void resetBaugrarteParameters(const double alpha, const double beta);

  // Resets the contact point by current kinematics of the robot. The kinematics
  // is passed through pinocchio::Data. Before calling this function, you have 
  // to update the kinematics (only related to configuraion) in pinocchio::Data.
  // Argments:
  //    data: The data including kinematics of the robot.
  void resetContactPointByCurrentKinematics(const pinocchio::Data& data);

  // Converts the double array of the contact force to joint forces.
  // Argments:
  //    contact_force: The contact force stored in double array.
  //    joint_force: The joint forces of the robot model. Resuts is stored in 
  //      this term.
  void computeJointForceFromContactForce(
      const Eigen::Vector3d& contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const;

  // Computes the 3xdimv contact Jacobian represented in the local coordinate 
  // of the contact frame. Before calling this function, you have to update the 
  // kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    J_contact: Contact jacobian is stored in this variable. size must be 
  //        3xdimv.
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          Eigen::MatrixXd& J_contact) const;

  // Computes the 3xdimv contact Jacobian represented in the local coordinate 
  // of the contact frame. Before calling this function, you have to update the 
  // kinematics of the model in pinocchio::Data.
  // The contact Jacobian is set in the block part of J_contacts with the 
  // indices of row and column begin at row_begin and column_begin.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    J_contacts: Contact jacobian is stored in this variable. size must be 
  //        at least 3xdimv.
  void getContactJacobian(const pinocchio::Model& model, pinocchio::Data& data, 
                          const unsigned int result_mat_row_begin, 
                          const unsigned int result_mat_col_begin,
                          Eigen::MatrixXd& J_contacts) const;

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_residual: The double array that the result is stored in. 
  void computeBaumgarteResidual(const pinocchio::Model& model, 
                                const pinocchio::Data& data, 
                                Eigen::Vector3d& baumgarte_residual) const;

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method. Before calling this function, you have 
  // to update the kinematics of the model in pinocchio::Data.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_residual: The double array that the result is stored in. 
  void computeBaumgarteResidual(const pinocchio::Model& model, 
                                const pinocchio::Data& data, 
                                const unsigned int result_vec_start_index,
                                Eigen::VectorXd& baumgarte_residual) const;

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. Before calling this 
  // function, you have to update the kinematics of the model in pinocchio::Data.
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
  // function, you have to update the kinematics of the model in pinocchio::Data.
  // The derivatives are set in the block part of baumgarte_partial_dq, 
  // baumgarte_partial_dv, and baumgarte_partial_da, with the indices of row 
  // and column begin at row_begin and column_begin.
  // Argments:
  //    model: pinocchio model of the robot.
  //    data: pinocchio data of the robot kinematics and dynamics.
  //    baumgarte_partial_dq: Partial of contact constraints with respect to 
  //      the generalized configuration. size must be at least 3xdimv.
  //    baumgarte_partial_dv: Partial of contact constraints with respect to 
  //      the generalized velocity. size must be at least 3xdimv.
  //    baumgarte_partial_da: Partial of contact constraints with respect to 
  //      the generalized acceleration. size must be at least 3xdimv.
  void computeBaumgarteDerivatives(const pinocchio::Model& model, 
                                   pinocchio::Data& data, 
                                   const unsigned int result_mat_row_begin, 
                                   const unsigned int result_mat_col_begin,
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
  unsigned int contact_frame_id() const;

  // Returns parent_joint_id, the index of the parent joint of the contact 
  // frame.
  unsigned int parent_joint_id() const;

  // Returns dimv, the dimensiton of the generalized velocity.
  unsigned int dimv() const;

  // Returns baumgarte_alpha, the weight parameter on the contact velocity in 
  // the Baumgarte's stabilization method.
  double baumgarte_alpha() const;

  // Returns baumgarte_alpha, the weight parameter on the contact position in 
  // the Baumgarte's stabilization method.
  double baumgarte_beta() const;
  
  // Returns the contact point.
  Eigen::Vector3d contact_point() const;

  // Returns SE(3) translation from the parent joint to the contact frame.
  pinocchio::SE3 jXf() const;

private:
  bool is_active_;
  unsigned int contact_frame_id_, parent_joint_id_, dimv_;
  double baumgarte_alpha_, baumgarte_beta_;
  Eigen::Vector3d contact_point_;
  pinocchio::SE3 jXf_;
  pinocchio::SE3::ActionMatrixType fXj_;
  Eigen::MatrixXd J_frame_, joint_v_partial_dq_, joint_a_partial_dq_, 
                  joint_a_partial_dv_, joint_a_partial_da_;
};

} // namespace idocp


#endif // IDOCP_POINT_CONTACT_HPP_ 