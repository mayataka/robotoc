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
  template <typename TangentVectorType, typename ConfigVectorType>
  void integrateConfiguration(
      const Eigen::MatrixBase<TangentVectorType>& v, 
      const double integration_length, 
      const Eigen::MatrixBase<ConfigVectorType>& q) const;

  template <typename ConfigVectorType1, typename TangentVectorType,  
            typename ConfigVectorType2>
  void integrateConfiguration(
      const Eigen::MatrixBase<ConfigVectorType1>& q, 
      const Eigen::MatrixBase<TangentVectorType>& v, 
      const double integration_length, 
      const Eigen::MatrixBase<ConfigVectorType2>& q_integrated) const;

  // Computes the difference of the two configurations at the its tangent 
  // velocity.
  // Argments: 
  //   q_plus: Configuration. Size must be dimq().
  //   q_minus: Configuration. Size must be dimq().
  //   difference: The resultant tangent vector of the configuration.
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename TangentVectorType>
  void subtractConfiguration(
      const Eigen::MatrixBase<ConfigVectorType1>& q_plus, 
      const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
      const Eigen::MatrixBase<TangentVectorType>& difference) const;

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
  template <typename ConfigVectorType, typename TangentVectorType,  
            typename MatrixType1, typename MatrixType2>
  void dIntegrateConfiguration(
      const Eigen::MatrixBase<ConfigVectorType>& q, 
      const Eigen::MatrixBase<TangentVectorType>& v, 
      const double integration_length, 
      const Eigen::MatrixBase<MatrixType1>& dIntegrate_dq,
      const Eigen::MatrixBase<MatrixType2>& dIntegrate_dv) const;

  // Differntiate the function of the subtraction of the configuration.
  // Argments: 
  //   q_plus: Configuration. Size must be dimq().
  //   q_minus: Configuration. Size must be dimq().
  //   dSubtract_dqplus: The partial derivative of the subtraction 
  //     q_plus - p_minus with respect to the q_plus.
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename MatrixType>
  void dSubtractdConfigurationPlus(
      const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
      const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
      const Eigen::MatrixBase<MatrixType>& dSubtract_dqplus) const;

  // Differntiate the function of the subtraction of the configuration.
  // Argments: 
  //   q_plus: Configuration. Size must be dimq().
  //   q_minus: Configuration. Size must be dimq().
  //   dSubtract_dqplus: The partial derivative of the subtraction 
  //     q_plus - p_minus with respect to the q_minus.
  template <typename ConfigVectorType1, typename ConfigVectorType2, 
            typename MatrixType>
  void dSubtractdConfigurationMinus(
      const Eigen::MatrixBase<ConfigVectorType1>& q_plus,
      const Eigen::MatrixBase<ConfigVectorType2>& q_minus,
      const Eigen::MatrixBase<MatrixType>& dSubtract_dminus) const;

  // Updates the kinematics of the robot. The frame placements, frame velocity,
  // frame acceleration, and the relevant Jacobians are calculated. After that, 
  // the each contact residual is updated.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2>
  void updateKinematics(const Eigen::MatrixBase<ConfigVectorType>& q, 
                        const Eigen::MatrixBase<TangentVectorType1>& v, 
                        const Eigen::MatrixBase<TangentVectorType2>& a);

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   residual: Vector where the result is stored. Size must be at least 3 and
  //     at most 3*max_point_contacts().
  template <typename VectorType>
  void computeBaumgarteResidual(
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   block_begin: The start index of the result.
  //   coeff: The coefficient of the result.
  //   residual: Vector where the result is stored. Size must be at least 3 and
  //     at most 3*max_point_contacts().
  template <typename VectorType>
  void computeBaumgarteResidual(
      const double coeff, 
      const Eigen::MatrixBase<VectorType>& baumgarte_residual) const;

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
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);

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
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void computeBaumgarteDerivatives(
      const double coeff, 
      const Eigen::MatrixBase<MatrixType1>& baumgarte_partial_dq, 
      const Eigen::MatrixBase<MatrixType2>& baumgarte_partial_dv, 
      const Eigen::MatrixBase<MatrixType3>& baumgarte_partial_da);
        
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
  template <typename VectorType>
  void setContactForces(const Eigen::MatrixBase<VectorType>& f);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3>
  void RNEA(const Eigen::MatrixBase<ConfigVectorType>& q, 
            const Eigen::MatrixBase<TangentVectorType1>& v, 
            const Eigen::MatrixBase<TangentVectorType2>& a, 
            const Eigen::MatrixBase<TangentVectorType3>& tau);

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
//   void RNEADerivatives(const Eigen::Ref<const Eigen::VectorXd>& q, 
//                        const Eigen::Ref<const Eigen::VectorXd>& v, 
//                        const Eigen::Ref<const Eigen::VectorXd>& a,
//                        Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_dq, 
//                        Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_dv, 
//                        Eigen::Ref<Eigen::MatrixXd> dRNEA_partial_da);
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename MatrixType1, 
            typename MatrixType2, typename MatrixType3>
  void RNEADerivatives(const Eigen::MatrixBase<ConfigVectorType>& q, 
                       const Eigen::MatrixBase<TangentVectorType1>& v, 
                       const Eigen::MatrixBase<TangentVectorType2>& a,
                       const Eigen::MatrixBase<MatrixType1>& dRNEA_partial_dq, 
                       const Eigen::MatrixBase<MatrixType2>& dRNEA_partial_dv, 
                       const Eigen::MatrixBase<MatrixType3>& dRNEA_partial_da);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to fext. This functin has to be called after calling 
  // updateKinematics().
  // Argments: 
  //   dRNEA_partial_dfext: The matrix where the result is stored. The size must 
  //      be at least 3*max_point_contacts() times dimv.
  template <typename MatrixType>
  void dRNEAPartialdFext(
      const Eigen::MatrixBase<MatrixType>& dRNEA_partial_dfext);

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
  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3,
            typename TangentVectorType4>
  void stateEquation(const Eigen::MatrixBase<ConfigVectorType>& q, 
                     const Eigen::MatrixBase<TangentVectorType1>& v, 
                     const Eigen::MatrixBase<TangentVectorType2>& tau, 
                     const Eigen::MatrixBase<TangentVectorType3>& dq,
                     const Eigen::MatrixBase<TangentVectorType4>& dv);

  // Generates feasible configuration randomly.
  // Argments:
  //   q: The generated configuration vector. Size must be dimq.  
  template <typename ConfigVectorType>
  void generateFeasibleConfiguration(
      const Eigen::MatrixBase<ConfigVectorType>& q) const;

  // Normalizes a configuration vector.
  // Argments:
  //   q: The normalized configuration vector. Size must be dimq.  
  template <typename ConfigVectorType>
  void normalizeConfiguration(
      const Eigen::MatrixBase<ConfigVectorType>& q) const;

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

  // Returns the dimensiton of joints.
  int dimJ() const;

  // Returns the maximum dimension of the contacts.
  int max_dimf() const;

  // Returns the dimension of the contacts.
  int dimf() const;

  // Returns the dimensiton of the generalized torques corresponding to the 
  // passive joints.
  int dim_passive() const;

  // Returns true if the robot has a floating base and false if not.
  bool has_floating_base() const;

  // Returns the maximum number of the contacts.
  int max_point_contacts() const;

  // Returns the number of the active point contacts.
  int num_active_point_contacts() const;

  // Returns true if contact[contact_index] is active. Returns false if 
  // contact[contact_index] is not active.
  bool is_contact_active(const int contact_index) const;

  // Prints the robot model.
  void printRobotModel() const;

private:

  void initializeJointLimits();

  pinocchio::Model model_;
  pinocchio::Data data_;
  std::string urdf_file_name_;
  std::vector<PointContact> point_contacts_;
  FloatingBase floating_base_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  int dimq_, dimv_, dimJ_, max_dimf_, dimf_, num_active_contacts_;
  std::vector<bool> is_each_contact_active_;
  Eigen::VectorXd joint_effort_limit_, joint_velocity_limit_,
                  lower_joint_position_limit_, upper_joint_position_limit_;
};

} // namespace idocp

#include "idocp/robot/robot.hxx"

#endif // IDOCP_ROBOT_HPP_ 