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
  ~StateEquation() = default;

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
  void correctCostateDirection(SplitDirection& d);

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

};

} // namespace robotoc 

#endif // ROBOTOC_STATE_EQUATION_HPP_