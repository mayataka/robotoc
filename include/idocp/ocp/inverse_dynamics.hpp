#ifndef IDOCP_INVERSE_DYNAMICS_HPP_
#define IDOCP_INVERSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"


namespace idocp {

class InverseDynamics {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  InverseDynamics(const Robot& robot) 
    : u_res_(Eigen::VectorXd::Zero(robot.dimv())),
      lu_(Eigen::VectorXd::Zero(robot.dimv())),
      lu_condensed_(Eigen::VectorXd::Zero(robot.dimv())),
      luu_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_dq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_dv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_da_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      du_df_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
      has_floating_base_(robot.has_floating_base()),
      dim_passive_(robot.dim_passive()),
      dimf_(robot.dimf()) {
  }

  InverseDynamics() 
    : u_res_(),
      lu_(),
      lu_condensed_(),
      luu_(),
      du_dq_(),
      du_dv_(),
      du_da_(),
      du_df_(),
      has_floating_base_(false),
      dim_passive_(0),
      dimf_(0) {
  }

  ~InverseDynamics() {
  }

  InverseDynamics(const InverseDynamics&) = default;

  InverseDynamics& operator=(const InverseDynamics&) = default;
 
  InverseDynamics(InverseDynamics&&) noexcept = default;

  InverseDynamics& operator=(InverseDynamics&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
  }

  inline void linearizeInverseDynamics(Robot& robot, const double dtau, 
                                       const SplitSolution& s) {
    if (dimf_ > 0) {
      robot.setContactForces(s.f);
    }
    robot.RNEA(s.q, s.v, s.a, u_res_);
    u_res_.noalias() -= s.u;
    robot.RNEADerivatives(s.q, s.v, s.a, du_dq_, du_dv_, du_da_);
    if (robot.dimf() > 0) {
      robot.dRNEAPartialdFext(du_df_);
    }
    lu_.noalias() -= dtau * s.beta; 
  }

  inline void linearizeFloatingBaseConstraint(const double dtau,
                                              const SplitSolution& s,
                                              KKTResidual& kkt_residual) {
    assert(dtau > 0);
    if (has_floating_base_) {
      lu_.head(dim_passive_).noalias() 
          += dtau * s.mu_active().tail(dim_passive_); 
      kkt_residual.C().tail(dim_passive_) = dtau * s.u.head(dim_passive_);
    }
  }

  inline void augmentInverseDynamicsDerivatives(
      const double dtau, const SplitSolution& s, 
      KKTResidual& kkt_residual) const {
    kkt_residual.lq().noalias() += dtau * du_dq_.transpose() * s.beta;
    kkt_residual.lv().noalias() += dtau * du_dv_.transpose() * s.beta;
    kkt_residual.la().noalias() += dtau * du_da_.transpose() * s.beta;
    if (dimf_ > 0) {
      kkt_residual.lf().noalias() += dtau * du_df_active_().transpose() * s.beta;
    }
  }

  inline void condenseFloatingBaseConstraint(const double dtau,
                                             KKTResidual& kkt_residual,
                                             KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    if (has_floating_base_) {
      kkt_residual.C().tail(dim_passive_).noalias()
          += dtau * u_res_.head(dim_passive_);
      kkt_matrix.Cq().bottomRows(dim_passive_) 
          = dtau * du_dq_.topRows(dim_passive_);
      kkt_matrix.Cv().bottomRows(dim_passive_) 
          = dtau * du_dv_.topRows(dim_passive_);
      kkt_matrix.Ca().bottomRows(dim_passive_) 
          = dtau * du_da_.topRows(dim_passive_);
      kkt_matrix.Cf().bottomRows(dim_passive_) 
          = dtau * du_df_active_().topRows(dim_passive_);
    }
  }

  inline void condenseInverseDynamics(KKTResidual& kkt_residual,
                                      KKTMatrix& kkt_matrix) {
    // condense Newton residual
    lu_condensed_ = lu_ + luu_ * u_res_;
    kkt_residual.lq().noalias() += du_dq_.transpose() * lu_condensed_;
    kkt_residual.lv().noalias() += du_dv_.transpose() * lu_condensed_;
    kkt_residual.la().noalias() += du_da_.transpose() * lu_condensed_;
    if (dimf_ > 0) {
      kkt_residual.lf().noalias() += du_df_active_().transpose() * lu_condensed_;
    }
    // condense Hessian
    kkt_matrix.Qaa() = du_da_.transpose() * luu_ * du_da_;
    if (dimf_ > 0) { 
      kkt_matrix.Qaf() = du_da_.transpose() * luu_ * du_df_active_(); 
    }
    kkt_matrix.Qaq() = du_da_.transpose() * luu_ * du_dq_;
    kkt_matrix.Qav() = du_da_.transpose() * luu_ * du_dv_;
    if (dimf_ > 0) {
      kkt_matrix.Qff() = du_df_active_().transpose() * luu_ * du_df_active_();
      kkt_matrix.Qfq() = du_df_active_().transpose() * luu_ * du_dq_;
      kkt_matrix.Qfv() = du_df_active_().transpose() * luu_ * du_dv_;
    }
    kkt_matrix.Qqq() = du_dq_.transpose() * luu_ * du_dq_;
    kkt_matrix.Qqv() = du_dq_.transpose() * luu_ * du_dv_;
    kkt_matrix.Qvv() = du_dv_.transpose() * luu_ * du_dv_;
  }

  inline void computeCondensedDirection(const double dtau, SplitDirection& d) {
    assert(dtau > 0);
    d.du = u_res_;
    d.du.noalias() += du_dq_ * d.dq();
    d.du.noalias() += du_dv_ * d.dv();
    d.du.noalias() += du_da_ * d.da();
    if (dimf_ > 0) {
      d.du.noalias() += du_df_active_() * d.df();
    }
    d.dbeta = (lu_  + luu_ * d.du) / dtau;
    if (has_floating_base_) {
      d.dbeta.head(dim_passive_).noalias() += d.dmu().tail(dim_passive_);
    }
  }

  inline double inverseDynamicsResidualL1Norm(const double dtau) const {
    return dtau * u_res_.lpNorm<1>();
  }

  inline double inverseDynamicsResidualL1Norm(Robot& robot, const double dtau,
                                              const SplitSolution& s) {
    if (dimf_ > 0) {
      robot.setContactForces(s.f);
    }
    robot.RNEA(s.q, s.v, s.a, u_res_);
    u_res_.noalias() -= s.u;
    return dtau * u_res_.lpNorm<1>();
  }

  inline double squaredKKTErrorNorm(const double dtau) const {
    return lu_.squaredNorm() + dtau * dtau * u_res_.squaredNorm();
  }

private:
  Eigen::VectorXd u_res_, lu_, lu_condensed_;
  Eigen::MatrixXd luu_, du_dq_, du_dv_, du_da_, du_df_;
  bool has_floating_base_;
  int dim_passive_, dimf_;

  inline Eigen::Ref<const Eigen::MatrixXd> du_df_active_() const {
    return du_df_.leftCols(dimf_);
  }

};

} // namespace idocp 


#endif // IDOCP_INVERSE_DYNAMICS_HPP_ 