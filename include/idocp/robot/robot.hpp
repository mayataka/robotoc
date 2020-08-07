#ifndef IDOCP_ROBOT_HPP_
#define IDOCP_ROBOT_HPP_

#include <string>
#include <map>
#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/aba.hpp"

#include "idocp/robot/point_contact.hpp"
#include "idocp/robot/floating_base.hpp"


namespace idocp {

class Robot {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Build the pinocchio model from urdf.
  Robot(const std::string& urdf_file_name);

  // Constructor. Build the pinocchio model from urdf. Build point contacts with
  // Baumgarte's stabilization parameters.
  Robot(const std::string& urdf_file_name, 
        const std::vector<int>& contact_frames, 
        const double baumgarte_weight_on_velocity, 
        const double baumgarte_weight_on_position);

  // Default constructor. 
  Robot();

  // Destructor. 
  ~Robot();

  // Use default copy constructor.
  Robot(const Robot&) = default;

  // Use default copy operator.
  Robot& operator=(const Robot&) = default;

  // Use default move constructor.
  Robot(Robot&&) noexcept = default;

  // Use default move assign operator.
  Robot& operator=(Robot&&) noexcept = default;

  // Build the pinocchio model from xml.
  void buildRobotModelFromXML(const std::string& xml);

  // Build the pinocchio model from xml. Build point contacts with Baumgarte's 
  // stabilization parameters.
  void buildRobotModelFromXML(const std::string& xml,
                              const std::vector<int>& contact_frames, 
                              const double baumgarte_weight_on_velocity, 
                              const double baumgarte_weight_on_position);

  // Integrates the generalized velocity, integration_length * v. 
  // The configuration q is then incremented.
  // Argments: 
  //   q: Configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  void integrateConfiguration(const Eigen::VectorXd& v, 
                              const double integration_length, 
                              Eigen::VectorXd& q) const;

  // Computes the difference of the two configurations at the its tangent 
  // velocity.
  // Argments: 
  //   q_plus: Configuration. Size must be dimq().
  //   q_minus: Configuration. Size must be dimq().
  //   difference: The resultant tangent vector of the configuration.
  void subtractConfiguration(const Eigen::VectorXd& q_plus, 
                             const Eigen::VectorXd& q_minus,
                             Eigen::VectorXd& difference) const;

  // Differntiate the function of integration the generalized velocity, 
  // integration_length * v, with respect to q and v.
  // Argments: 
  //   q: Configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  //   dIntegrate_dq: The partial derivative of the integration with respect to
  //     the configuratioh q.
  //   dIntegrate_dv: The partial derivative of the integration with respect to
  //     the generalized velocity v.
  void dIntegrateConfiguration(const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v,
                               const double integration_length,
                               Eigen::MatrixXd& dIntegrate_dq,
                               Eigen::MatrixXd& dIntegrate_dv) const;

  // Computes the Jacobian of the configuration computed at its tangent 
  // velocity.
  // Argments: 
  //   q: Configuration. Size must be dimq().
  void computeConfigurationJacobian(const Eigen::VectorXd& q);

  // Transforms the gradient evaluated at the configuration space into that 
  // evaluated at its tangent space. Before calling this function, 
  // computeConfigurationJacobian() must be called. 
  // Argments: 
  //   gradient_at_configuration: The gradient evaluated at the configuration 
  //     space. Size must be dimq().
  //   gradient_at_configuration: The gradient evaluated at the tangent
  //     space. Size must be dimv().
  void computeTangentGradient(const Eigen::VectorXd& gradient_at_configuration, 
                              Eigen::VectorXd& gradient_at_tangent) const;

  // Transforms the Hessian evaluated at the configuration space into that 
  // evaluated at its tangent space. Before calling this function,
  // computeConfigurationJacobian() must be called. 
  // Argments: 
  //   hessian_at_configuration: The Hessian evaluated at the configuration 
  //     space. Size must be dimq() x dimq().
  //   hessian_at_tangent: The Hessian evaluated at the tangent
  //     space. Size must be dimv() x dimv().
  void computeTangentHessian(const Eigen::MatrixXd& hessian_at_configuration, 
                             Eigen::MatrixXd& hessian_at_tangent) const;

  // Transforms the Hessian evaluated at the configuration space into that 
  // evaluated at its tangent space and add into another Hessian. Before 
  // calling this function, computeConfigurationJacobian() must be called. 
  // Argments: 
  //   hessian_at_configuration: The Hessian evaluated at the configuration 
  //     space. Size must be dimq() x dimq().
  //   coeff: The coefficient at augmenting the hessian.
  //   augmented_hessian_at_tangent: The augmented Hessian the transformed
  //    hessian_at_configuration is added. Size must be dimv() x dimv().
  void augmentTangentHessian(const Eigen::MatrixXd& hessian_at_configuration, 
                             const double coeff,
                             Eigen::MatrixXd& augmented_hessian_at_tangent) const;

  // Updates the kinematics of the robot. The frame placements, frame velocity,
  // frame acceleration, and the relevant Jacobians are calculated. After that, 
  // the each contact residual is updated.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  void updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Eigen::VectorXd& a);

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   residual: Vector where the result is stored. Size must be at least 3 and
  //     at most 3*max_point_contacts().
  void computeBaumgarteResidual(const int block_begin, 
                                Eigen::VectorXd& baumgarte_residual);

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   block_begin: The start index of the result.
  //   coeff: The coefficient of the result.
  //   residual: Vector where the result is stored. Size must be at least 3 and
  //     at most 3*max_point_contacts().
  void computeBaumgarteResidual(const int block_begin, const double coeff, 
                                Eigen::VectorXd& baumgarte_residual);

  // Computes the product of a vector and the derivatives of the contact 
  // constriants represented by Baumgarte's stabilization method. 
  // Before calling this function, updateKinematics() must be called. 
  // Argments: 
  //   dBaumgarte_partial_dq: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  //   dBaumgarte_partial_dv: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  //   dBaumgarte_partial_da: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  void computeBaumgarteDerivatives(const int block_rows_begin,
                                   Eigen::MatrixXd& dBaumgarte_partial_dq, 
                                   Eigen::MatrixXd& dBaumgarte_partial_dv,
                                   Eigen::MatrixXd& dBaumgarte_partial_da);

  // Computes the product of a vector and the derivatives of the contact 
  // constriants represented by Baumgarte's stabilization method. 
  // Before calling this function, updateKinematics() must be called. 
  // Argments: 
  //   block_rows_begin: The start index of the block rows where result stored.
  //   coeff: The coefficient of the result.
  //   dBaumgarte_partial_dq: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  //   dBaumgarte_partial_dv: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  //   dBaumgarte_partial_da: The matrix where the result is stored. The number 
  //     of columns must be dimv. The number of rows must be at least 3 and 
  //     at most 3*max_point_contacts().
  void computeBaumgarteDerivatives(const int block_rows_begin, 
                                   const double coeff,
                                   Eigen::MatrixXd& dBaumgarte_partial_dq, 
                                   Eigen::MatrixXd& dBaumgarte_partial_dv,
                                   Eigen::MatrixXd& dBaumgarte_partial_da);

  // Sets the contact points.
  void setContactPoints(const std::vector<Eigen::Vector3d>& contact_points);

  // Sets all the contact points by using current kinematics.
  void setContactPointsByCurrentKinematics();

  // Activate and deactivate the each contact.
  //   is_each_contact_active: containts the bool variables representing 
  //     wheather each contact is active or not.
  void setContactStatus(const std::vector<bool>& is_each_contact_active);

  // Set the stack of the contact forces. Before calling this function, call
  // setContactStatus().
  //   fext: The stack of the contact forces represented in the local frame.
  //      The size must be at most max_dimf().
  void setContactForces(const Eigen::VectorXd& fext);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const Eigen::VectorXd& a, Eigen::VectorXd& tau);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, and a.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   dRNEA_partial_dq: The matrix where the result is stored. The size must 
  //      be dimv times dimv.
  //   dRNEA_partial_dv: The matrix where the result is stored. The size must 
  //      be dimv times dimv.
  //   dRNEA_partial_da: The matrix where the result is stored. The size must 
  //      be dimv times dimv.
  void RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, 
                       Eigen::MatrixXd& dRNEA_partial_dq, 
                       Eigen::MatrixXd& dRNEA_partial_dv, 
                       Eigen::MatrixXd& dRNEA_partial_da);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to fext. This functin has to be called after calling 
  // updateKinematics().
  // Argments: 
  //   dRNEA_partial_dfext: The matrix where the result is stored. The size must 
  //      be at least 3*max_point_contacts() times dimv.
  void dRNEAPartialdFext(Eigen::MatrixXd& dRNEA_partial_dfext);

  // void impulse(const Eigen::VectorXd& q, const Eigen::VectorXd& v_before,
  //              const Eigen::MatrixXd& J_contacts, Eigen::VectorXd& v_after);

  // Computes the state equation, i.e., the corresponding generalized velocity 
  // and the generalized acceleration under given q, v, tau.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  //   dq: Returned generalized velocity. Size must be dimv.
  //   dv: Returned generalized acceleration. Size must be dimv.
  void stateEquation(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& tau, Eigen::VectorXd& dq, 
                     Eigen::VectorXd& dv);

  // Substitutes zero in the generalized torques tau corresponding to the 
  // passive joints.
  // Argments:
  //   torques: The generalized torques for fully actuated system. The size 
  //     must be dimv.
  void setPassiveTorques(Eigen::VectorXd& torques) const;

  // Generates feasible configuration randomly.
  // Argments:
  //   q: The generated configuration vector. Size must be dimq.  
  void generateFeasibleConfiguration(Eigen::VectorXd& q) const;

  // Normalizes a configuration vector.
  // Argments:
  //   q: The normalized configuration vector. Size must be dimq.  
  void normalizeConfiguration(Eigen::VectorXd& q) const;

  // Returns the effort limit of each joints.
  Eigen::VectorXd jointEffortLimit() const;

  // Returns the joint velocity limit of each joints.
  Eigen::VectorXd jointVelocityLimit() const;

  // Returns the lower limit of the position of each joints.
  Eigen::VectorXd lowerJointPositionLimit() const;

  // Returns the upper limit of the position of each joints.
  Eigen::VectorXd upperJointPositionLimit() const;

  // Sets the effort limit of each joints.
  void setJointEffortLimit(const Eigen::VectorXd& joint_effort_limit);

  // Sets the joint velocity limit of each joints.
  void setJointVelocityLimit(const Eigen::VectorXd& joint_velocity_limit);

  // Sets the lower limit of the position of each joints.
  void setLowerJointPositionLimit(
      const Eigen::VectorXd& lower_joint_position_limit);

  // Sets the upper limit of the position of each joints.
  void setUpperJointPositionLimit(
      const Eigen::VectorXd& upper_joint_position_limit);

  // Returns the dimensiton of the generalized configuration.
  int dimq() const;

  // Returns the dimensiton of the generalized velocity.
  int dimv() const;

  // Returns the dimension of the contacts.
  int dimf() const;

  // Returns the maximum dimension of the contacts.
  int max_dimf() const;

  // Returns true if the robot has a floating base and false if not.
  bool has_floating_base() const;

  // Returns the dimensiton of the generalized torques corresponding to the 
  // passive joints.
  int dim_passive() const;

  std::vector<int> passive_joint_indices() const;

  // Returns the maximum number of the contacts.
  int max_point_contacts() const;

  // Returns true if contact[contact_index] is active. Returns false if 
  // contact[contact_index] is not active.
  bool is_contact_active(const int contact_index) const;

  // Prints the robot model.
  void printRobotModel() const;

private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  std::string urdf_file_name_;
  std::vector<PointContact> point_contacts_;
  FloatingBase floating_base_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  int dimq_, dimv_, dimf_, max_dimf_;
  std::vector<bool> is_each_contact_active_;
  Eigen::VectorXd joint_effort_limit_, joint_velocity_limit_,
                  lower_joint_position_limit_, upper_joint_position_limit_;
  Eigen::MatrixXd configuration_jacobian_;
};

} // namespace idocp


#endif // IDOCP_ROBOT_HPP_ 