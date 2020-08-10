#ifndef IDOCP_PARNMPC_LINEARIZER_HPP_
#define IDOCP_PARNMPC_LINEARIZER_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints_interface.hpp"


namespace idocp {
namespace parnmpclinearizer {

inline void linearizeStageCost(Robot& robot, 
                               std::shared_ptr<CostFunctionInterface>& cost, 
                               CostFunctionData& cost_data,
                               const double t, const double dtau, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, 
                               const Eigen::VectorXd& a, 
                               const Eigen::VectorXd& f, 
                               const Eigen::VectorXd& u, 
                               Eigen::Ref<Eigen::VectorXd> lq, 
                               Eigen::Ref<Eigen::VectorXd> lv, 
                               Eigen::Ref<Eigen::VectorXd> la, 
                               Eigen::Ref<Eigen::VectorXd> lf, 
                               Eigen::Ref<Eigen::VectorXd> lu) {
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
    robot.computeConfigurationJacobian(q, cost_data.configuration_jacobian);
  }
  cost->lq(robot, cost_data, t, dtau, q, v, a, lq);
  cost->lv(robot, cost_data, t, dtau, q, v, a, lv);
  cost->la(robot, cost_data, t, dtau, q, v, a, la);
  cost->lf(robot, cost_data, t, dtau, f, lf);
  cost->lu(robot, cost_data, t, dtau, u, lu);
}


inline void linearizeDynamics(Robot& robot, const double dtau,
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const Eigen::VectorXd& a, 
                              const Eigen::VectorXd& f, 
                              const Eigen::VectorXd& u, 
                              const Eigen::VectorXd& q_prev, 
                              const Eigen::VectorXd& v_prev, 
                              Eigen::VectorXd& kkt_res, Eigen::VectorXd& u_res, 
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
  robot.subtractConfiguration(q_prev, q, kkt_res.head(robot.dimv()));
  kkt_res.head(robot.dimv()).noalias() += dtau * v;
  kkt_res.segment(robot.dimv(), robot.dimv()) = v_prev - v + dtau * a;
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
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const Eigen::VectorXd& a, 
                                 const Eigen::VectorXd& u, 
                                 const Eigen::VectorXd& u_res, 
                                 const Eigen::MatrixXd& du_dq, 
                                 const Eigen::MatrixXd& du_dv, 
                                 const Eigen::MatrixXd& du_da, 
                                 const Eigen::MatrixXd& du_df, 
                                 Eigen::VectorXd& kkt_res, 
                                 Eigen::MatrixXd& kkt_mat, 
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
  const int dimf = robot.dimf();
  const int dimv = robot.dimv();
  const int constraint_block_rows_begin = 2*dimv;
  const int constraint_block_cols_begin = 2*dimv + dimf + dim_passive;
  if (dimf > 0) {
    robot.computeBaumgarteResidual(constraint_block_rows_begin, dtau, kkt_res);
    robot.computeBaumgarteDerivatives(dim_passive, dtau, Cq, Cv, Ca);
    kkt_mat.block(constraint_block_rows_begin, constraint_block_cols_begin, 
                  dimf, dimv) = Ca.topRows(dimf);
    kkt_mat.block(constraint_block_rows_begin, 
                  constraint_block_cols_begin+dimv+dimf, 
                  dimf, dimv) = Cq.topRows(dimf);
    kkt_mat.block(constraint_block_rows_begin, 
                  constraint_block_cols_begin+2*dimv+dimf, 
                  dimf, dimv) = Cv.topRows(dimf);
  }
  if (dim_passive > 0) {
    kkt_res.segment(constraint_block_rows_begin+dimf, dim_passive) 
        = dtau * (u.head(dim_passive)+u_res.head(dim_passive));
    kkt_mat.block(constraint_block_rows_begin+dimf, constraint_block_cols_begin, 
                  dim_passive, dimv) = dtau * du_dq.topRows(dim_passive);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+dimv, 
                  dim_passive, dimf) = dtau * du_df.topLeftCorner(dim_passive, 
                                                                  dimf);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+dimv+dimf, 
                  dim_passive, dimv) = dtau * du_dq.topLeftCorner(dim_passive, 
                                                                  dimv);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+2*dimv+dimf, 
                  dim_passive, dimv) = dtau * du_dv.topLeftCorner(dim_passive, 
                                                                  dimv);
  }
}


inline void factorizeKKTResidual(const Eigen::VectorXd& q_res, 
                                 const Eigen::VectorXd& v_res,
                                 const Eigen::VectorXd& C_res,
                                 const Eigen::VectorXd& lq, 
                                 const Eigen::VectorXd& lv, 
                                 const Eigen::VectorXd& la, 
                                 const Eigen::VectorXd& lf, 
                                 Eigen::VectorXd& kkt_res) {
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dimc = robot.dimf() + robot.dim_passive();
  kkt_res.segment(0, dimv) = q_res;
  kkt_res.segment(dimv, dimv) = v_res;
  kkt_res.segment(2*dimv, dimc) = C_res;
  kkt_res.segment(2*dimv+dimc, dimv) = la;
  kkt_res.segment(3*dimv+dimc, dimf) = lf.head(dimf);
  kkt_res.segment(3*dimv+dimc+dimf, dimv) = lq;
  kkt_res.segment(4*dimv+dimc+dimf, dimv) = lv;
}


inline void factorizeKKTMatrix(const Eigen::MatrixXd& lqq, 
                               const Eigen::MatrixXd& lvv, 
                               const Eigen::MatrixXd& laa, 
                               const Eigen::MatrixXd& lff, 
                               Eigen::MatrixXd& kkt_mat) {
  const int constraint_block_rows_begin = 2*dimv;
  const int constraint_block_cols_begin = 2*dimv + dimf + dim_passive;
  if (dimf > 0) {
    robot.computeBaumgarteResidual(constraint_block_rows_begin, dtau, kkt_res);
    robot.computeBaumgarteDerivatives(dim_passive, dtau, Cq, Cv, Ca);
    kkt_mat.block(constraint_block_rows_begin, constraint_block_cols_begin, 
                  dimf, dimv) = Ca.topRows(dimf);
    kkt_mat.block(constraint_block_rows_begin, 
                  constraint_block_cols_begin+dimv+dimf, 
                  dimf, dimv) = Cq.topRows(dimf);
    kkt_mat.block(constraint_block_rows_begin, 
                  constraint_block_cols_begin+2*dimv+dimf, 
                  dimf, dimv) = Cv.topRows(dimf);
  }
  if (dim_passive > 0) {
    kkt_res.segment(constraint_block_rows_begin+dimf, dim_passive) 
        = dtau * (u.head(dim_passive)+u_res.head(dim_passive));
    kkt_mat.block(constraint_block_rows_begin+dimf, constraint_block_cols_begin, 
                  dim_passive, dimv) = dtau * du_dq.topRows(dim_passive);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+dimv, 
                  dim_passive, dimf) = dtau * du_df.topLeftCorner(dim_passive, 
                                                                  dimf);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+dimv+dimf, 
                  dim_passive, dimv) = dtau * du_dq.topLeftCorner(dim_passive, 
                                                                  dimv);
    kkt_mat.block(constraint_block_rows_begin+dimf, 
                  constraint_block_cols_begin+2*dimv+dimf, 
                  dim_passive, dimv) = dtau * du_dv.topLeftCorner(dim_passive, 
                                                                  dimv);
  }
}


} // namespace parnmpclinearizer
} // namespace idocp 


#endif // IDOCP_PARNMPC_LINEARIZER_HPP_