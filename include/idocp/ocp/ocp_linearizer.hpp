#ifndef IDOCP_OCP_LINEARIZER_HPP_
#define IDOCP_OCP_LINEARIZER_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints_interface.hpp"


namespace idocp {
namespace ocplinearizer {

inline void linearizeStageCost(Robot& robot, 
                              std::shared_ptr<CostFunctionInterface>& cost, 
                              CostFunctionData& cost_data,
                              const double t, const double dtau, 
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const Eigen::VectorXd& a, 
                              const Eigen::VectorXd& u, 
                              const Eigen::VectorXd& f, Eigen::VectorXd& lq, 
                              Eigen::VectorXd& lv, Eigen::VectorXd& la, 
                              Eigen::VectorXd& lu, Eigen::VectorXd& lf) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  assert(f.size() == robot.max_dimf());
  assert(lq.size() == robot.dimv());
  assert(lv.size() == robot.dimv());
  assert(la.size() == robot.dimv());
  assert(lu.size() == robot.dimv());
  assert(lf.size() == robot.max_dimf());
  if (robot.has_floating_base()) {
    robot.computeConfigurationJacobian(q);
  }
  cost->lq(robot, cost_data, t, dtau, q, v, a, lq);
  cost->lv(robot, cost_data, t, dtau, q, v, a, lv);
  cost->la(robot, cost_data, t, dtau, q, v, a, la);
  cost->lu(robot, cost_data, t, dtau, u, lu);
  cost->lf(robot, cost_data, t, dtau, f, lf);
}


inline void linearizeDynamics(Robot& robot, const double dtau,
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const Eigen::VectorXd& a, 
                              const Eigen::VectorXd& u, 
                              const Eigen::VectorXd& f, 
                              const Eigen::VectorXd& q_next, 
                              const Eigen::VectorXd& v_next, 
                              Eigen::VectorXd& q_res,
                              Eigen::VectorXd& v_res, Eigen::VectorXd& u_res, 
                              Eigen::MatrixXd& du_dq, Eigen::MatrixXd& du_dv, 
                              Eigen::MatrixXd& du_da, Eigen::MatrixXd& du_df) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  assert(f.size() == robot.max_dimf());
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  assert(q_res.size() == robot.dimv());
  assert(v_res.size() == robot.dimv());
  assert(u_res.size() == robot.dimv());
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(du_df.rows() == robot.dimv());
  assert(du_df.cols() == robot.max_dimf());
  robot.subtractConfiguration(q, q_next, q_res);
  q_res.noalias() += dtau * v;
  v_res = v + dtau * a - v_next;
  robot.setContactForces(f);
  robot.RNEA(q, v, a, u_res);
  u_res.noalias() -= u;
  robot.RNEADerivatives(q, v, a, du_dq, du_dv, du_da);
  robot.dRNEAPartialdFext(du_df);
}


inline void linearizeConstraints(Robot& robot, const double dtau,
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const Eigen::VectorXd& a, 
                                 const Eigen::VectorXd& u, 
                                 const Eigen::VectorXd& u_res, 
                                 const Eigen::MatrixXd& du_dq, 
                                 const Eigen::MatrixXd& du_dv, 
                                 const Eigen::MatrixXd& du_da, 
                                 const Eigen::MatrixXd& du_df, 
                                 Eigen::VectorXd& C_res, 
                                 Eigen::MatrixXd& Cq, Eigen::MatrixXd& Cv, 
                                 Eigen::MatrixXd& Ca, Eigen::MatrixXd& Cf) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
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
  assert(C_res.size() == robot.dim_passive()+robot.max_dimf());
  assert(Cq.rows() == robot.dim_passive()+robot.max_dimf());
  assert(Cq.cols() == robot.dimv());
  assert(Cv.rows() == robot.dim_passive()+robot.max_dimf());
  assert(Cv.cols() == robot.dimv());
  assert(Ca.rows() == robot.dim_passive()+robot.max_dimf());
  assert(Ca.cols() == robot.dimv());
  assert(Cf.rows() == robot.dim_passive());
  assert(Cf.cols() == robot.max_dimf());
  const int dim_passive = robot.dim_passive();
  C_res.head(dim_passive) = dtau * (u.head(dim_passive)+u_res.head(dim_passive));
  Cq.topRows(dim_passive) = dtau * du_dq.topRows(dim_passive);
  Cv.topRows(dim_passive) = dtau * du_dv.topRows(dim_passive);
  Ca.topRows(dim_passive) = dtau * du_da.topRows(dim_passive);
  const int dimf = robot.dimf();
  Cf.leftCols(dimf) = dtau * du_df.topLeftCorner(dim_passive, dimf);
  robot.computeBaumgarteResidual(dim_passive, dtau, C_res);
  robot.computeBaumgarteDerivatives(dim_passive, dtau, Cq, Cv, Ca);
}

} // namespace ocplinearizer
} // namespace idocp 


#endif // IDOCP_OCP_LINEARIZER_HPP_