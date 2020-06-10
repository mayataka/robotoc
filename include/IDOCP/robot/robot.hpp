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

#include "robot/point_contact.hpp"
#include "robot/passive_joints.hpp"


namespace idocp {

class Robot {
public:
  Robot(const std::string& urdf_file_name);

  Robot(const std::string& urdf_file_name, 
        const std::vector<unsigned int>& contact_frames, 
        const double baumgarte_weight_on_position, 
        const double baumgarte_weight_on_velocity);

  Robot();

  // Destructor. 
  ~Robot();

  // Use default copy constructor.
  Robot(const Robot& other) = default;

  // Use default copy operator.
  Robot& operator=(const Robot& other) = default;

  // Integrates the generalized velocity, integration_length * v. 
  // The generalized configuration q is then incremented.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  void integrateConfiguration(const Eigen::VectorXd& v, 
                              const double integration_length, 
                              Eigen::VectorXd& q) const;

  // Integrates the generalized velocity, integration_length * v. 
  // The generalized configuration q is then incremented.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  void differenceConfiguration(const Eigen::VectorXd& q_plus, 
                               const Eigen::VectorXd& q_minus,
                               Eigen::VectorXd& difference) const;

  void dIntegrateConfiguration(const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v,
                               const double integration_length,
                               Eigen::MatrixXd& dIntegrate_dq,
                               Eigen::MatrixXd& dIntegrate_dv) const;

  // Updates the kinematics of the robot. The frame placements, frame velocity,
  // frame acceleration, and the relevant Jacobians are calculated. After that, 
  // the each contact residual is updated.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  void updateKinematics(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Eigen::VectorXd& a);

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   residual: The array where the result is stored. The must be dimf.
  void computeBaumgarteResidual(Eigen::VectorXd& baumgarte_residual) const;

  // Computes the product of a vector and the derivatives of the contact 
  // constriants represented by Baumgarte's stabilization method. 
  // Before calling this function, updateKinematics() must be called.
  // Argments: 
  //   vec: A vector multiplied to partial derivatives of the contact 
  //      contstraints represented by Baumgarte's method. 
  //      The size is assumed to be dimf.
  //   dBaum_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dBaum_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dBaum_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void computeBaumgarteDerivatives(Eigen::MatrixXd& dBaumgarte_partial_dq, 
                                   Eigen::MatrixXd& dBaumgarte_partial_dv,
                                   Eigen::MatrixXd& dBaumgarte_partial_da);
                                      

  // Sets the stack of the generalized forces represented in the world frame.
  // The array is converted into joint forces.
  //   fext: The stack of the contact forces represented in the world frame.
  //      The size is assumed to be dimf.
  void setActiveContacts(const std::vector<bool>& is_each_contact_active, 
                         const Eigen::VectorXd& fext);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const Eigen::VectorXd& a, Eigen::VectorXd& tau);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, and a, and returns products of them and vec.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, 
                       Eigen::MatrixXd& dRNEA_partial_dq, 
                       Eigen::MatrixXd& dRNEA_partial_dv, 
                       Eigen::MatrixXd& dRNEA_partial_da);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, and a, and returns products of them and vec.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void dRNEAPartialdFext(Eigen::MatrixXd& dRNEA_partial_dfext);

  // Substitutes zero in the generalized torques tau corresponding to the 
  // passive joints.
  // Argments:
  //   tau: The generalized torque for fully actuated system. The size is dimv.
  void setPassiveTorques(Eigen::VectorXd& tau) const;

  // Calculates the violation of torques corresponding to the passive joints 
  // under given generalized torques tau.
  // Argments:
  //   tau: The generalized torque for fully actuated system. The size is dimv.
  //   violation: The residual of the constraints of the zero torques. The size
  //      is dim_passive.
  void passiveConstraintViolation(const Eigen::VectorXd& tau, 
                                  Eigen::VectorXd& violation) const;

  Eigen::VectorXd jointEffortLimit() const;

  Eigen::VectorXd jointVelocityLimit() const;

  Eigen::VectorXd lowerJointPositionLimit() const;

  Eigen::VectorXd upperJointPositionLimit() const;

  // Returns the dimensiton of the generalized configuration.
  unsigned int dimq() const;

  // Returns the dimensiton of the generalized velocity.
  unsigned int dimv() const;

  // Returns the dimensiton of the generalized torques corresponding to the 
  // passive joints.
  unsigned int dim_passive() const;

  // Returns the maximum number of the point contacts.
  unsigned int max_point_contacts() const;

private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  std::string urdf_file_name_;
  std::vector<PointContact> point_contacts_;
  PassiveJoints passive_joints_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  unsigned int dimq_, dimv_;
};

} // namespace idocp


#endif // IDOCP_OCP_ROBOT_HPP_ 