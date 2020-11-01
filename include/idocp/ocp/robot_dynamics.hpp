#ifndef IDOCP_ROBOT_DYNAMICS_HPP_
#define IDOCP_ROBOT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/contact_dynamics.hpp"
#include "idocp/ocp/non_contact_dynamics.hpp"


namespace idocp {

class RobotDynamics {
public:
  RobotDynamics(const Robot& robot);

  RobotDynamics();

  ~RobotDynamics();

  RobotDynamics(const RobotDynamics&) = default;

  RobotDynamics& operator=(const RobotDynamics&) = default;
 
  RobotDynamics(RobotDynamics&&) noexcept = default;

  RobotDynamics& operator=(RobotDynamics&&) noexcept = default;

  void linearizeRobotDynamics(Robot& robot, const ContactStatus& contact_status, 
                              const double dtau, const SplitSolution& s, 
                              KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  void condenseRobotDynamics(Robot& robot, const ContactStatus& contact_status,
                             const double dtau, const SplitSolution& s, 
                             KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const double dtau, 
                                       const KKTMatrix& kkt_matrix, 
                                       const KKTResidual& kkt_residual, 
                                       SplitDirection& d);

  void computeCondensedDualDirection(const double dtau, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual, 
                                     SplitDirection& d);

  void computeRobotDynamicsResidual(Robot& robot, 
                                    const ContactStatus& contact_status,
                                    const double dtau, const SplitSolution& s, 
                                    KKTResidual& kkt_residual);

  double l1NormRobotDynamicsResidual(const double dtau) const;

  double squaredNormRobotDynamicsResidual(const double dtau) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_, has_active_contacts_;
  ContactDynamics contact_dynamics_;
  NonContactDynamics non_contact_dynamics_;

};

} // namespace idocp 

#include "idocp/ocp/robot_dynamics.hxx"

#endif // IDOCP_ROBOT_DYNAMICS_HPP_