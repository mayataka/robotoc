#ifndef IDOCP_STATE_EQUATION_HPP_
#define IDOCP_STATE_EQUATION_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class StateEquation {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  StateEquation(const Robot& robot) 
    : dsubtract_dqminus_(),
      dsubtract_dqplus_() {
    if (robot.has_floating_base()) {
      dsubtract_dqminus_.resize(robot.dimv(), robot.dimv());
      dsubtract_dqplus_.resize(robot.dimv(), robot.dimv());
      dsubtract_dqminus_.setZero();
      dsubtract_dqplus_.setZero();
    }
  }

  StateEquation() 
    : dsubtract_dqminus_(),
      dsubtract_dqplus_() {
  }

  ~StateEquation() {
  }

  StateEquation(const StateEquation&) = default;

  StateEquation& operator=(const StateEquation&) = default;
 
  StateEquation(StateEquation&&) noexcept = default;

  StateEquation& operator=(StateEquation&&) noexcept = default;

  inline void linearizeStateAndCostateEquation(
      Robot& robot, const double dtau, const Eigen::VectorXd& q_prev, 
      const Eigen::VectorXd& v_prev, const SplitSolution& s, 
      const Eigen::VectorXd& lmd_next, const Eigen::VectorXd& gmm_next, 
      const Eigen::VectorXd& q_next, KKTResidual& kkt_residual, 
      KKTMatrix& kkt_matrix) {
    assert(dtau > 0);
    assert(q_prev.size() == robot.dimq());
    assert(v_prev.size() == robot.dimv());
    assert(lmd_next.size() == robot.dimv());
    assert(gmm_next.size() == robot.dimv());
    assert(q_next.size() == robot.dimq());
    robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
    kkt_residual.Fq().noalias() += dtau * s.v;
    kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq());
      robot.dSubtractdConfigurationPlus(s.q, q_next, kkt_residual.dsubtract_dq);
      kkt_residual.lq().noalias() 
          += kkt_matrix.Fqq().transpose() * s.lmd 
              + kkt_residual.dsubtract_dq.transpose() * lmd_next;
    }
    else {
      kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
      kkt_residual.lq().noalias() += lmd_next - s.lmd;
    }
    kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm + gmm_next;
    kkt_residual.la().noalias() += dtau * s.gmm;
    kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
    kkt_matrix.Fvv() = - Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
    kkt_matrix.Fva() = dtau * Eigen::MatrixXd::Identity(robot.dimv(), 
                                                        robot.dimv());
  }

  inline void linearizeStateAndCostateEquationTerminal(
      Robot& robot, const double dtau, const Eigen::VectorXd& q_prev, 
      const Eigen::VectorXd& v_prev, const SplitSolution& s, 
      KKTResidual& kkt_residual, KKTMatrix& kkt_matrix) {
    assert(dtau > 0);
    assert(q_prev.size() == robot.dimq());
    assert(v_prev.size() == robot.dimv());
    robot.subtractConfiguration(q_prev, s.q, kkt_residual.Fq());
    kkt_residual.Fq().noalias() += dtau * s.v;
    kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq());
      kkt_residual.lq().noalias() 
          += kkt_matrix.Fqq().transpose() * s.lmd;
    }
    else {
      kkt_residual.lq().noalias() += - s.lmd;
    }
    kkt_residual.lv().noalias() += dtau * s.lmd - s.gmm;
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

private:
  Eigen::MatrixXd dsubtract_dqminus_, dsubtract_dqplus_;

};

} // namespace idocp 


#endif // IDOCP_STATE_EQUATION_HPP_