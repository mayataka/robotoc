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

class ParNMPCLinearizer
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  ParNMPCLinearizer(const Robot& robot) 
    : dsubtract_dqminus_(),
      dsubtract_dqplus_(),
      phiq_(Eigen::VectorXd::Zero(robot.dimv())),
      phiv_(Eigen::VectorXd::Zero(robot.dimv())) {
    if (robot.has_floating_base()) {
      dsubtract_dqminus_.resize(robot.dimv(), robot.dimv());
      dsubtract_dqplus_.resize(robot.dimv(), robot.dimv());
      dsubtract_dqminus_.setZero();
      dsubtract_dqplus_.setZero();
    }
  }

  ParNMPCLinearizer() 
    : dsubtract_dqminus_(),
      dsubtract_dqplus_(),
      phiq_(),
      phiv_() {
  }

  ~ParNMPCLinearizer() {
  }

  ParNMPCLinearizer(const ParNMPCLinearizer&) = default;

  ParNMPCLinearizer& operator=(const ParNMPCLinearizer&) = default;
 
  ParNMPCLinearizer(ParNMPCLinearizer&&) noexcept = default;

  ParNMPCLinearizer& operator=(ParNMPCLinearizer&&) noexcept = default;

  inline void linearizeStageCost(Robot& robot, 
                                 const std::shared_ptr<CostFunction>& cost, 
                                 CostFunctionData& cost_data, const double t, 
                                 const double dtau, const SplitSolution& s,
                                 KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    cost->lq(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.lq());
    cost->lv(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.lv());
    cost->la(robot, cost_data, t, dtau, s.q, s.v, s.a, kkt_residual.la());
    cost->lf(robot, cost_data, t, dtau, s.f, kkt_residual.lf());
  }

  inline void linearizeStateEquation(
      Robot& robot, const double dtau, 
      const Eigen::Ref<const Eigen::VectorXd>& q_prev, 
      const Eigen::Ref<const Eigen::VectorXd>& v_prev, const SplitSolution& s, 
      const Eigen::Ref<const Eigen::VectorXd>& lmd_next, 
      const Eigen::Ref<const Eigen::VectorXd>& gmm_next, 
      const Eigen::Ref<const Eigen::VectorXd>& q_next, 
      KKTResidual& kkt_residual, KKTMatrix& kkt_matrix) {
    assert(dtau > 0);
    assert(q_prev.size() == robot.dimq());
    assert(v_prev.size() == robot.dimv());
    robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
    kkt_residual.Fq().noalias() += dtau * s.v;
    kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus_);
      robot.dSubtractdConfigurationPlus(s.q, q_next, dsubtract_dqplus_);
      kkt_residual.lq().noalias() 
          += dsubtract_dqminus_.transpose() * s.lmd 
              + dsubtract_dqplus_.transpose() * lmd_next;
    }
    else {
      kkt_residual.lq().noalias() += lmd_next - s.lmd;
    }
    kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm + gmm_next;
    kkt_residual.la().noalias() += dtau * s.gmm;
    if (robot.has_floating_base()) {
      kkt_matrix.Fqq() = dsubtract_dqminus_;
    }
    else {
      kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
    }
    kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
    kkt_matrix.Fvv() = - Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
    kkt_matrix.Fva() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
  }

  inline void linearizeTerminalCost(Robot& robot, 
                                    const std::shared_ptr<CostFunction>& cost, 
                                    CostFunctionData& cost_data, const double t, 
                                    const SplitSolution& s) {
    cost->phiq(robot, cost_data, t, s.q, s.v, phiq_);
    cost->phiv(robot, cost_data, t, s.q, s.v, phiv_);
  }

  inline void linearizeStateEquation(
      Robot& robot, const double dtau,
      const Eigen::Ref<const Eigen::VectorXd>& q_prev, 
      const Eigen::Ref<const Eigen::VectorXd>& v_prev, const SplitSolution& s, 
      KKTResidual& kkt_residual, KKTMatrix& kkt_matrix) {
    assert(dtau > 0);
    assert(q_prev.size() == robot.dimq());
    assert(v_prev.size() == robot.dimv());
    robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
    kkt_residual.Fq().noalias() += dtau * s.v;
    kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus_);
      kkt_residual.lq().noalias() 
          += dsubtract_dqminus_.transpose() * s.lmd + phiq_;
    }
    else {
      kkt_residual.lq().noalias() += phiq_ - s.lmd;
    }
    kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm + phiv_;
    kkt_residual.la().noalias() += dtau * s.gmm;
    if (robot.has_floating_base()) {
      kkt_matrix.Fqq() = dsubtract_dqminus_;
    }
    else {
      kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
    }
    kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
    kkt_matrix.Fvv() = - Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
    kkt_matrix.Fva() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
  }

  inline void linearizeContactConstraints(Robot& robot, const double dtau,
                                          KKTResidual& kkt_residual,
                                          KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
    robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                      kkt_matrix.Ca());
  }

private:
  Eigen::MatrixXd dsubtract_dqminus_, dsubtract_dqplus_;
  Eigen::VectorXd phiq_, phiv_;

};

} // namespace idocp 


#endif // IDOCP_PARNMPC_LINEARIZER_HPP_