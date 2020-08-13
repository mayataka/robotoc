#ifndef IDOCP_PARNMPC_LINEARIZER_HPP_
#define IDOCP_PARNMPC_LINEARIZER_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {
namespace parnmpclinearizer {

inline void linearizeStageCost(Robot& robot, std::shared_ptr<CostFunction>& cost, 
                               CostFunctionData& cost_data, const double t, 
                               const double dtau, const SplitSolution& s,
                               KKTResidual& kkt_residual,
                               Eigen::Ref<Eigen::VectorXd> lu) {
  assert(dtau > 0);
  assert(lu.size() == robot.dimv());
  cost->lq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.lq());
  cost->lv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.lv());
  cost->la(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.la());
  cost->lf(robot, cost_data, t, dtau, s.f, kkt_residual.lf());
  cost->lu(robot, cost_data, t, dtau, s.u, lu);
}


inline void linearizeDynamics(Robot& robot, const double dtau,
                              const Eigen::Ref<const Eigen::VectorXd>& q_prev, 
                              const Eigen::Ref<const Eigen::VectorXd>& v_prev, 
                              const SplitSolution& s, KKTResidual& kkt_residual,
                              Eigen::Ref<Eigen::VectorXd> u_res,
                              Eigen::Ref<Eigen::MatrixXd> du_dq,
                              Eigen::Ref<Eigen::MatrixXd> du_dv,
                              Eigen::Ref<Eigen::MatrixXd> du_da,
                              Eigen::Ref<Eigen::MatrixXd> du_df) {
  assert(dtau > 0);
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
  robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
  kkt_residual.Fq().noalias() += dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  if (robot.dimf() > 0) {
    robot.setContactForces(s.f);
  }
  robot.RNEA(s.q, s.v, s.a, u_res);
  u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  if (robot.dimf() > 0) {
    robot.dRNEAPartialdFext(du_df);
  }
}


inline void linearizeConstraints(Robot& robot, const double dtau,
                                 const SplitSolution& s,
                                 const Eigen::Ref<const Eigen::VectorXd>& u_res, 
                                 const Eigen::Ref<const Eigen::MatrixXd>& du_dq, 
                                 const Eigen::Ref<const Eigen::MatrixXd>& du_dv, 
                                 const Eigen::Ref<const Eigen::MatrixXd>& du_da, 
                                 const Eigen::Ref<const Eigen::MatrixXd>& du_df, 
                                 KKTResidual& kkt_residual,
                                 KKTMatrix& kkt_matrix) {
  assert(dtau > 0);
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
    kkt_residual.C().tail(robot.dim_passive()) = dtau * (s.u.head(robot.dim_passive())+u_res.head(robot.dim_passive()));
    kkt_matrix.Cq().bottomRows(robot.dim_passive()) = dtau * du_dq.topRows(robot.dim_passive());
    kkt_matrix.Cv().bottomRows(robot.dim_passive()) = dtau * du_dv.topRows(robot.dim_passive());
    kkt_matrix.Ca().bottomRows(robot.dim_passive()) = dtau * du_da.topRows(robot.dim_passive());
    kkt_matrix.Cf().bottomRows(robot.dim_passive()) = dtau * du_df.topLeftCorner(robot.dim_passive(), robot.dimf());
  }
}

} // namespace parnmpclinearizer
} // namespace idocp 


#endif // IDOCP_PARNMPC_LINEARIZER_HPP_