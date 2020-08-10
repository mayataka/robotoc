#ifndef IDOCP_PARNMPC_LINEARIZER_HPP_
#define IDOCP_PARNMPC_LINEARIZER_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {
namespace parnmpclinearizer {

inline void linearizeStageCost(Robot& robot, 
                               std::shared_ptr<CostFunctionInterface>& cost, 
                               CostFunctionData& cost_data,
                               const double t, const double dtau, 
                               const Eigen::Ref<const Eigen::VectorXd> q, 
                               const Eigen::Ref<const Eigen::VectorXd> v, 
                               const Eigen::Ref<const Eigen::VectorXd> a, 
                               const Eigen::Ref<const Eigen::VectorXd> f, 
                               const Eigen::Ref<const Eigen::VectorXd> u, 
                               KKTResidual& kkt_residual,
                               Eigen::Ref<Eigen::VectorXd> Qu) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(f.size() == robot.max_dimf());
  assert(u.size() == robot.dimv());
  if (robot.has_floating_base()) {
    robot.computeConfigurationJacobian(q, cost_data.configuration_jacobian);
  }
  cost->lq(robot, cost_data, t, dtau, q, v, a, kkt_residual.Qq());
  cost->lv(robot, cost_data, t, dtau, q, v, a, kkt_residual.Qv());
  cost->la(robot, cost_data, t, dtau, q, v, a, kkt_residual.Qa());
  cost->lf(robot, cost_data, t, dtau, f, kkt_residual.Qf());
  cost->lu(robot, cost_data, t, dtau, u, Qu);
}


inline void linearizeDynamics(Robot& robot, const double dtau,
                              const Eigen::Ref<const Eigen::VectorXd> q, 
                              const Eigen::Ref<const Eigen::VectorXd> v, 
                              const Eigen::Ref<const Eigen::VectorXd> a, 
                              const Eigen::Ref<const Eigen::VectorXd> f, 
                              const Eigen::Ref<const Eigen::VectorXd> u, 
                              const Eigen::Ref<const Eigen::VectorXd> q_prev, 
                              const Eigen::Ref<const Eigen::VectorXd> v_prev, 
                              KKTResidual& kkt_residual,
                              Eigen::Ref<Eigen::VectorXd> u_res,
                              Eigen::Ref<Eigen::MatrixXd> du_dq,
                              Eigen::Ref<Eigen::MatrixXd> du_dv,
                              Eigen::Ref<Eigen::MatrixXd> du_da,
                              Eigen::Ref<Eigen::MatrixXd> du_df) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(f.size() == robot.max_dimf());
  assert(u.size() == robot.dimv());
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  assert(u_res.size() == robot.dimv());
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(du_df.rows() == robot.dimv());
  assert(du_df.cols() == robot.max_dimf());
  robot.subtractConfiguration(q_prev, q, kkt_residual.Fq());
  kkt_residual.Fq().head(robot.dimv()).noalias() += dtau * v;
  kkt_residual.Fv() = v_prev - v + dtau * a;
  if (robot.dimf() > 0) {
    robot.setContactForces(f);
  }
  robot.RNEA(q, v, a, u_res);
  u_res.noalias() -= u;
  robot.RNEADerivatives(q, v, a, du_dq, du_dv, du_da);
  if (robot.dimf() > 0) {
    robot.dRNEAPartialdFext(du_df);
  }
}


inline void linearizeConstraints(Robot& robot, const double dtau,
                                 const Eigen::Ref<const Eigen::VectorXd> u, 
                                 const Eigen::Ref<const Eigen::VectorXd> u_res, 
                                 const Eigen::Ref<const Eigen::MatrixXd> du_dq, 
                                 const Eigen::Ref<const Eigen::MatrixXd> du_dv, 
                                 const Eigen::Ref<const Eigen::MatrixXd> du_da, 
                                 const Eigen::Ref<const Eigen::MatrixXd> du_df, 
                                 KKTResidual& kkt_residual,
                                 KKTMatrix& kkt_matrix) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(u_res.size() == robot.dimv());
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(du_df.rows() == robot.dimv());
  assert(du_df.cols() == robot.max_dimf());
  if (robot.dimf() > 0) {
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
    robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                      kkt_matrix.Ca());
  }
  if (robot.dim_passive() > 0) {
    kkt_residual.C().tail(robot.dim_passive())
        = dtau * (u.head(robot.dim_passive())+u_res.head(robot.dim_passive()));
    kkt_matrix.Cq().bottomRows(robot.dim_passive()) 
        = dtau * du_dq.topRows(robot.dim_passive());
    kkt_matrix.Cv().bottomRows(robot.dim_passive()) 
        = dtau * du_dv.topRows(robot.dim_passive());
    kkt_matrix.Ca().bottomRows(robot.dim_passive()) 
        = dtau * du_da.topRows(robot.dim_passive());
    kkt_matrix.Cf().bottomRows(robot.dim_passive()) 
        = dtau * du_df.topRows(robot.dim_passive());
  }
}

} // namespace parnmpclinearizer
} // namespace idocp 


#endif // IDOCP_PARNMPC_LINEARIZER_HPP_