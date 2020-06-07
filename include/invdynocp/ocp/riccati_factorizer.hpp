#ifndef INVDYNOCP_RICCATI_FACTORIZER_HPP_
#define INVDYNOCP_RICCATI_FACTORIZER_HPP_

#include <vector>

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace invdynocp {

class RiccatiFactorizer {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  RiccatiFactorizer(const Robot& robot, const CostFunctionInterface* cost,
                    const ConstraintsInterface* constraints);

  ~RiccatiFactorizer();
 
  // Copy constructor.
  RiccatiFactorizer(const RiccatiFactorizer& other) = default;

  // Copy operator.
  RiccatiFactorizer& operator=(const RiccatiFactorizer& other) = default;

  // without contacts 
  void linearizeOCP(Robot& robot, const double dtau, const Eigen::VectorXd& x, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& mu,
                    const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& x_next);

  void linearizeOCP(Robot& robot, const double dtau, const Eigen::VectorXd& x, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& fext, 
                    const Eigen::VectorXd& mu, const Eigen::VectorXd& nu,
                    const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& x_next);

  void factorizeRiccatiMatrix(Eigen::MatrixXd& A, Eigen::MatrixXd& B);

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

private:
  unsigned int dim_Tx_, dim_afext_;
  Eigen::VectorXd u_, x_res_, lmd_res_, afext_res_, C_res_;
  Eigen::MatrixXd Fx_, Fafext_, Hx_, Hafext_;
  Eigen::MatrixXd du_dq_, du_dv_, du_dafext_;

};

} // namespace invdynocp


#endif // INVDYNOCP_RICCATI_FACTORIZER_HPP_