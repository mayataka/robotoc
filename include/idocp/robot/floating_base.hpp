#ifndef IDOCP_FLOATING_BASE_HPP_
#define IDOCP_FLOATING_BASE_HPP_

#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"


namespace idocp {

class FloatingBase {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized, e.g., by pinocchio::buildModel() or 
  //      pinocchio::buildModelFromXML().
  FloatingBase(const pinocchio::Model& model);

  // Default constructor. 
  FloatingBase();
 
  // Destructor. 
  ~FloatingBase();

  // Use dafault copy constructor.
  FloatingBase(const FloatingBase&) = default;

  // Use dafule copy operator.
  FloatingBase& operator=(const FloatingBase&) = default;

  // Substitutes zero in the generalized torques tau corresponding to the 
  // passive joints.
  // Argments:
  //   torques: The generalized torques for fully actuated system. The size 
  //     must be dimv.
  void setPassiveTorques(Eigen::VectorXd& torques) const;

  // Calculates the violation of torques corresponding to the passive joints 
  // under given generalized torques.
  // Argments:
  //   torques: The generalized torque for fully actuated system. The size must   
  //     be dimv.
  //   violation: The residual of the constraints of the zero torques. The size
  //      must be dim_passive.
  void computePassiveConstraintViolation(const Eigen::VectorXd& torques, 
                                         Eigen::VectorXd& violation) const;

  // Returns the dimension of the torques correspoinding to the passive joints.
  int dim_passive() const;

  // Returns the indices of the virtual passive joints corresponding to the 
  // floating base.
  std::vector<int> passive_joint_indices() const;

  // Returns true if the robot has a floating base and false if not.
  bool has_floating_base() const;

private:
  bool has_floating_base_;
  std::vector<int> passive_joint_indices_;
  int dimv_;
};

} // namespace idocp


#endif // IDOCP_FLOATING_BASE_HPP_