#ifndef IDOCP_INVERSE_DYNAMICS_HPP_
#define IDOCP_INVERSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class InverseDynamics {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  InverseDynamics(const Robot& robot);

  InverseDynamics();

  ~InverseDynamics();

  InverseDynamics(const InverseDynamics&) = default;

  InverseDynamics& operator=(const InverseDynamics&) = default;
 
  InverseDynamics(InverseDynamics&&) noexcept = default;

  InverseDynamics& operator=(InverseDynamics&&) noexcept = default;

  void linearizeInverseDynamics(Robot& robot, const double dtau, 
                                const SplitSolution& s,
                                KKTResidual& kkt_residual);

  void condenseInverseDynamics(KKTMatrix& kkt_matrix, 
                               KKTResidual& kkt_residual);

  void condenseEqualityConstraint(const double dtau, KKTMatrix& kkt_matrix,
                                  KKTResidual& kkt_residual) const;

  void computeCondensedDirection(const double dtau, 
                                 const KKTMatrix& kkt_matrix, 
                                 const KKTResidual& kkt_residual, 
                                 SplitDirection& d);

  double violationL1Norm(const double dtau, 
                         const KKTResidual& kkt_residual) const;

  double violationL1Norm(Robot& robot, const double dtau, 
                         const SplitSolution& s, 
                         KKTResidual& kkt_residual) const;

private:
  Eigen::VectorXd lu_condensed_;
  Eigen::MatrixXd du_dq_, du_dv_, du_da_, du_df_;
  bool has_floating_base_, has_active_contacts_;
  int dim_passive_, dimf_;

  void setContactStatus(const Robot& robot);

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  du_df_active_();

  const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  du_df_active_() const;

};

} // namespace idocp 

#include "idocp/ocp/inverse_dynamics.hxx"

#endif // IDOCP_INVERSE_DYNAMICS_HPP_ 