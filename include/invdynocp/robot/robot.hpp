#ifndef INVDYN_OCP_ROBOT_HPP_
#define INVDYN_OCP_ROBOT_HPP_

#include <string>
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


namespace invdynocp {

class Robot {
public:
  // Constructor. Creates workspace for pinocchio.
  Robot(const std::string& urdf_file_name, 
        const unsigned int max_point_contacts);

  // Destructor. 
  ~Robot();

  // Copy constructor.
  Robot(const Robot& other);

  // Copy operator.
  Robot& operator=(const Robot& other);

  // Move constructor.
  Robot(const Robot&& other) noexcept;

  // Copy operator.
  Robot& operator=(const Robot&& other) noexcept;

  // Integrates the generalized velocity, integration_length * v. 
  // The generalized configuration q is then incremented.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  void integrateConfiguration(const Eigen::VectorXd& v, 
                              const double integration_length,
                              Eigen::VectorXd& q);

  // Integrates the generalized velocity, integration_length * v. 
  // The initial configuration is evolves to q_plus.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq().
  //   v: Generalized velocity. Size must be dimv().
  //   integration_length: The length of the integration.
  //   q_plus: Generalized configuration after integration. Size must be dimq().
  void integrateConfiguration(const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const double integration_length, 
                              Eigen::VectorXd& q_plus);

  // Computes difference of the two configurations represented in the tangent 
  // space. The tangent vector from q_minus to q_plus is stored in difference.
  // Argments: 
  //   q_plus: Generalized configuration. Size must be dimq().
  //   q_minus: Generalized configuration. Size must be dimq().
  //   difference: Generalized velocity representing the difference. Size must 
  //      be dimv().
  void differenceConfigurations(const Eigen::VectorXd& q_plus, 
                                const Eigen::VectorXd& q_minus, 
                                Eigen::VectorXd& difference);

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
  void computeBaumgarteResidual(Eigen::VectorXd& baumgarte_residual);

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
  void setFext(const Eigen::VectorXd& fext);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const Eigen::VectorXd& a, Eigen::VectorXd& tau);

  // Computes generalized torques tau corresponding to given q, v, a, and fext.
  // External forces are set by setFext() before calling this function.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   fext: The stack of the contact forces. Size must be dimf.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const Eigen::VectorXd& a, const Eigen::VectorXd& fext, 
            Eigen::VectorXd& tau);

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
  void RNEADerivatives(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, const Eigen::VectorXd& fext, 
                       Eigen::MatrixXd& dRNEA_partial_dq, 
                       Eigen::MatrixXd& dRNEA_partial_dv, 
                       Eigen::MatrixXd& dRNEA_partial_da_and_fext);

  // Adds the point contact to the robot.
  // If there is a contact that have contact_frame_id, the contact does not 
  //    increased. Otherwise, a contact that has contact_frame_id is added.
  // Argments:
  //    contact_frame_id: The frame index of the contact. 
  //    baumgarte_alpha: The weight parameter of the Baumgrate's stabilization
  //      method
  //    baumgarte_beta: The weight parameter of the Baumgrate's stabilization
  //      method
  void addPointContact(const unsigned int contact_frame_id, 
                       const double baumgarte_alpha, 
                       const double baumgarte_beta);

  // Removes the point contact from the robot. If there is no contact that has 
  //    contact_frame_id, this function does not do anything.
  // Argments:
  //    contact_frame_id: The frame index of the contact. 
  void removePointContact(const unsigned int contact_frame_id);

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

  // Calculates the product of a vector and the derivative of the residual of 
  // the constrains on the passive joints. In usually, the vector must be 
  // the corresponding Lagrange multiplier.
  // Argments:
  //   vec: A vector whose size is dim_passive.
  //   added_vec: A vector which vec is added. The size must be Robot::dimv().
  void passiveConstraintsDerivative(Eigen::MatrixXd& derivative) const;

  // Returns the dimensiton of the generalized configuration.
  unsigned int dimq() const;

  // Returns the dimensiton of the generalized velocity.
  unsigned int dimv() const;

  // Returns the dimensiton of the contacts that is equal to the dimension of 
  //    the contact forces.
  unsigned int dimf() const;

  // Returns the maximum dimensiton of the contacts that is equal to the 
  // dimension of the contact forces.
  unsigned int dimfmax() const;

  // Returns the dimensiton of the generalized torques corresponding to the 
  // passive joints.
  unsigned int dim_passive() const;

  // Returns the maximum number of the point contacts.
  unsigned int max_point_contacts() const;

  // Returns the name of the URDF file name.
  std::string urdf_file_name() const;

private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  const std::string urdf_file_name_;
  std::vector<PointContact> point_contacts_;
  PassiveJoints passive_joints_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  unsigned int dimq_, dimv_, max_point_contacts_;

};

} // namespace invdynocp


#endif // INVDYN_OCP_ROBOT_HPP_