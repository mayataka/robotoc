#ifndef INVDYNOCP_RICCATI_HPP_
#define INVDYNOCP_RICCATI_HPP_

#include <vector>

#include "Eigen/Core"


namespace invdynocp {

class Riccati {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  Riccati(const unsigned int dimx, const unsigned int dimu, 
          const unsigned dimc);

  ~Riccati();
 
  // Destructor. 
  ~PassiveJoints();

  void setPassiveTorques(Eigen::VectorXd& tau) const;

  void computePassiveConstraintViolation(const Eigen::VectirXd& torques, 
                                         Eigen::VectirXd& violation) const;

  void computePassiveConstraintDerivative(Eigen::MatrixXd& derivative) const;

  // Returns the dimension of the torques correspoinding to the passive joints.
  int dim_passive() const;

  // Prohibits copy constructor.
  PassiveJoints(const PassiveJoints&) = default;

  // Prohibits copy operator.
  PassiveJoints& operator=(const PassiveJoints&) = default;

private:
  unsigned int dimx_, dimu_, dimc_;

};

} // namespace invdynocp


#endif // INVDYNOCP_RICCATI_HPP_