#ifndef INVDYNOCP_PASSIVE_JOINTS_HPP_
#define INVDYNOCP_PASSIVE_JOINTS_HPP_

#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"


namespace invdynocp {

class PassiveJoints {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  PassiveJoints(const pinocchio::Model& model);

  PassiveJoints();
 
  // Destructor. 
  ~PassiveJoints();

  // Use dafault copy constructor.
  PassiveJoints(const PassiveJoints&) = default;

  // Use dafule copy operator.
  PassiveJoints& operator=(const PassiveJoints&) = default;


  void setPassiveTorques(Eigen::VectorXd& tau) const;

  void computePassiveConstraintViolation(const Eigen::VectorXd& torques, 
                                         Eigen::VectorXd& violation) const;

  void computePassiveConstraintDerivative(Eigen::MatrixXd& derivative) const;

  // Returns the dimension of the torques correspoinding to the passive joints.
  unsigned int dim_passive() const;

  std::vector<unsigned int> passive_torque_indices() const;

private:
  std::vector<unsigned int> passive_torque_indices_;
};

} // namespace invdynocp


#endif // INVDYNOCP_PASSIVE_JOINTS_HPP_