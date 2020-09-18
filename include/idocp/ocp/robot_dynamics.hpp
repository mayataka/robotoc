#ifndef IDOCP_ROBOT_DYNAMICS_HPP_
#define IDOCP_ROBOT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


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

  void setContactStatus(const Robot& robot);

  void augmentRobotDynamics(Robot& robot, const double dtau, 
                            const SplitSolution& s, KKTMatrix& kkt_matrix, 
                            KKTResidual& kkt_residual);

  void condenseRobotDynamics(Robot& robot, const double dtau,
                             const SplitSolution& s, KKTMatrix& kkt_matrix, 
                             KKTResidual& kkt_residual);

  void computeCondensedDirection(const double dtau, 
                                 const KKTMatrix& kkt_matrix, 
                                 const KKTResidual& kkt_residual, 
                                 SplitDirection& d);

  double violationL1Norm(Robot& robot, const double dtau, 
                         const SplitSolution& s, 
                         KKTResidual& kkt_residual) const;

  static double computeViolationL1Norm(Robot& robot, const double dtau, 
                                       const SplitSolution& s, 
                                       KKTResidual& kkt_residual);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4, typename MatrixType5, typename MatrixType6>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& da_dq,
                            const Eigen::MatrixBase<MatrixType2>& da_dv,
                            const Eigen::MatrixBase<MatrixType3>& df_dq,
                            const Eigen::MatrixBase<MatrixType4>& df_dv,
                            const Eigen::MatrixBase<MatrixType5>& Kuq,
                            const Eigen::MatrixBase<MatrixType6>& Kuv) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd lu_condensed_;
  Eigen::MatrixXd du_dq_, du_dv_, du_da_, du_df_full_, 
                  Quu_du_dq_, Quu_du_dv_, Quu_du_da_, Quu_du_df_full_;
  bool has_floating_base_, has_active_contacts_;
  int dimf_;
  static constexpr int kDimFloatingBase = 6;

  void linearizeInverseDynamics(Robot& robot, const SplitSolution& s, 
                                KKTResidual& kkt_residual);

  static void linearizeContactConstraint(Robot& robot, const double dtau,
                                         KKTMatrix& kkt_matrix, 
                                         KKTResidual& kkt_residual);

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  du_df_();

  const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  du_df_() const;

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  Quu_du_df_();

  const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  Quu_du_df_() const;

};

} // namespace idocp 

#include "idocp/ocp/robot_dynamics.hxx"

#endif // IDOCP_ROBOT_DYNAMICS_HPP_