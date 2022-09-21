#ifndef ROBOTOC_STATE_EQUATION_HPP_
#define ROBOTOC_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class StateEquation
/// @brief State equation. Only represent kinematic relation between stages.
///
class StateEquation {
public:
  ///
  /// @brief Constructs a state equation.
  /// @param[in] robot Robot model. 
  ///
  StateEquation(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  StateEquation();
  
  ///
  /// @brief Destructor. 
  ///
  ~StateEquation();

  ///
  /// @brief Default copy constructor. 
  ///
  StateEquation(const StateEquation&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  StateEquation& operator=(const StateEquation&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateEquation(StateEquation&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateEquation& operator=(StateEquation&&) noexcept = default;

  ///
  /// @brief Computes the residual in the state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step. 
  /// @param[in] s Solution at the current time stage. 
  /// @param[in] q_next Configuration at the next time stage. 
  /// @param[in] v_next Generalized velocity at the next time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  static void evalStateEquation(const Robot& robot, const double dt, 
                                const SplitSolution& s, 
                                const Eigen::VectorXd& q_next, 
                                const Eigen::VectorXd& v_next, 
                                SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  static void linearizeStateEquation(const Robot& robot, const double dt, 
                                     const Eigen::VectorXd& q_prev, 
                                     const SplitSolution& s, 
                                     const SplitSolution& s_next, 
                                     SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the linearized state equation using the Jacobian of the 
  /// Lie group. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step. 
  /// @param[in] s Solution at the current stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix at the current time stage. 
  /// @param[in, out] kkt_residual Split KKT residual at the current time stage. 
  ///
  void correctLinearizedStateEquation(const Robot& robot, const double dt, 
                                      const SplitSolution& s, 
                                      const SplitSolution& s_next, 
                                      SplitKKTMatrix& kkt_matrix, 
                                      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the costate direction using the Jacobian of the Lie group. 
  /// @param[in, out] d Split direction. 
  ///
  void correctCostateDirection(SplitDirection& d)  {
    if (has_floating_base_) {
      Fq_tmp_ = Fqq_prev_inv_.transpose() * d.dlmdgmm.template head<6>();
      d.dlmdgmm.template head<6>() = - Fq_tmp_;
    }
  }

  ///
  /// @brief Computes the initial state direction using the result of  
  /// StateEquation::linearizeStateEquationAlongLieGroup().
  /// @param[in] robot Robot model. 
  /// @param[in] q0 Initial configuration. 
  /// @param[in] v0 Initial generalized velocity. 
  /// @param[in] s0 Split solution at the initial stage. 
  /// @param[in, out] d0 Split direction at the initial stage. 
  ///
  void computeInitialStateDirection(const Robot& robot, 
                                    const Eigen::VectorXd& q0, 
                                    const Eigen::VectorXd& v0, 
                                    const SplitSolution& s0, 
                                    SplitDirection& d0) const;

private:
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_, Fqq_tmp_;  
  Eigen::VectorXd Fq_tmp_;
  SE3JacobianInverse se3_jac_inverse_;
  bool has_floating_base_;

  template <typename SplitSolutionType>
  static void linearizeStateEquation_impl(const Robot& robot, const double dt, 
                                          const Eigen::VectorXd& q_prev, 
                                          const SplitSolution& s, 
                                          const SplitSolutionType& s_next, 
                                          SplitKKTMatrix& kkt_matrix, 
                                          SplitKKTResidual& kkt_residual) {
    assert(dt > 0);
    assert(q_prev.size() == robot.dimq());
    evalStateEquation(robot, dt, s, s_next.q, s_next.v, kkt_residual);
    if (robot.hasFloatingBase()) {
      robot.dSubtractConfiguration_dqf(s.q, s_next.q, kkt_matrix.Fqq());
      robot.dSubtractConfiguration_dq0(q_prev, s.q, kkt_matrix.Fqq_prev);
      kkt_residual.lq().template head<6>().noalias() 
          += kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
                * s_next.lmd.template head<6>();
      kkt_residual.lq().template head<6>().noalias() 
          += kkt_matrix.Fqq_prev.template topLeftCorner<6, 6>().transpose() 
                * s.lmd.template head<6>();
      kkt_residual.lq().tail(robot.dimv()-6).noalias() 
          += s_next.lmd.tail(robot.dimv()-6) - s.lmd.tail(robot.dimv()-6);
    }
    else {
      kkt_matrix.Fqq().diagonal().fill(1.);
      kkt_residual.lq().noalias() += s_next.lmd - s.lmd;
    }
    kkt_matrix.Fqv().diagonal().fill(dt);
    kkt_residual.lv().noalias() += dt * s_next.lmd + s_next.gmm - s.gmm;
    kkt_residual.la.noalias() += dt * s_next.gmm;
    // STO sensitivities
    kkt_residual.h += s_next.lmd.dot(s.v);
    kkt_residual.h += s_next.gmm.dot(s.a);
    kkt_matrix.hv().noalias() += s_next.lmd;
    kkt_matrix.ha.noalias()   += s_next.gmm;
    kkt_matrix.fq() = s.v;
    kkt_matrix.fv() = s.a;
  }

  template <typename SplitSolutionType>
  void correctLinearizedStateEquation_impl(const Robot& robot, const double dt, 
                                           const SplitSolution& s, 
                                           const SplitSolutionType& s_next, 
                                           SplitKKTMatrix& kkt_matrix, 
                                           SplitKKTResidual& kkt_residual) {
    if (has_floating_base_) {
      assert(dt > 0);
      se3_jac_inverse_.compute(kkt_matrix.Fqq_prev, Fqq_prev_inv_);
      robot.dSubtractConfiguration_dq0(s.q, s_next.q, kkt_matrix.Fqq_prev);
      se3_jac_inverse_.compute(kkt_matrix.Fqq_prev, Fqq_inv_);
      Fqq_tmp_ = kkt_matrix.Fqq().template topLeftCorner<6, 6>();
      kkt_matrix.Fqq().template topLeftCorner<6, 6>().noalias() = - Fqq_inv_ * Fqq_tmp_;
      kkt_matrix.Fqv().template topLeftCorner<6, 6>() = - dt * Fqq_inv_;
      Fq_tmp_  = kkt_residual.Fq().template head<6>();
      kkt_residual.Fq().template head<6>().noalias() = - Fqq_inv_ * Fq_tmp_;
      Fq_tmp_ = kkt_matrix.fq().template head<6>();
      kkt_matrix.fq().template head<6>().noalias() = - Fqq_inv_ * Fq_tmp_;
    }
  }

};

} // namespace robotoc 

#endif // ROBOTOC_STATE_EQUATION_HPP_