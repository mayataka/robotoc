#ifndef ROBOTOC_IMPULSE_STATE_EQUATION_HPP_
#define ROBOTOC_IMPULSE_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3_jacobian_inverse.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class ImpulseStateEquation
/// @brief State equation. Only represent kinematic relation
/// between stages.
///
class ImpulseStateEquation {
public:
  ///
  /// @brief Constructs an impulse state equation.
  /// @param[in] robot Robot model. 
  ///
  ImpulseStateEquation(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  ImpulseStateEquation();
  
  ///
  /// @brief Destructor. 
  ///
  ~ImpulseStateEquation() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseStateEquation(const ImpulseStateEquation&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseStateEquation& operator=(const ImpulseStateEquation&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseStateEquation(ImpulseStateEquation&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseStateEquation& operator=(ImpulseStateEquation&&) noexcept = default;

  ///
  /// @brief Computes the residual in the impulse state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] q_next Configuration at the next time stage. 
  /// @param[in] v_next Generalized velocity at the next time stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  static void evalStateEquation(const Robot& robot, 
                                const SplitSolution& s, 
                                const Eigen::VectorXd& q_next, 
                                const Eigen::VectorXd& v_next, 
                                SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the impulse state equation. 
  /// @param[in] robot Robot model. 
  /// @param[in] q_prev Configuration at the previous time stage. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
  /// stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  static void linearizeStateEquation(const Robot& robot, 
                                     const Eigen::VectorXd& q_prev, 
                                     const SplitSolution& s, 
                                     const SplitSolution& s_next, 
                                     SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the linearized state equation using the Jacobian of the 
  /// Lie group. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Solution at the current impulse stage. 
  /// @param[in] s_next Solution at the next time stage. 
  /// @param[in, out] kkt_matrix Impulse split KKT matrix at the current impulse 
  /// stage. 
  /// @param[in, out] kkt_residual Impulse split KKT residual at the current 
  /// impulse stage. 
  ///
  void correctLinearizedStateEquation(const Robot& robot, 
                                      const SplitSolution& s, 
                                      const SplitSolution& s_next, 
                                      SplitKKTMatrix& kkt_matrix, 
                                      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Corrects the costate direction using the Jacobian of the Lie group. 
  /// @param[in, out] d Split direction. 
  ///
  void correctCostateDirection(SplitDirection& d);

private:
  Eigen::MatrixXd Fqq_inv_, Fqq_prev_inv_, Fqq_tmp_;  
  Eigen::VectorXd Fq_tmp_;
  SE3JacobianInverse se3_jac_inverse_;
  bool has_floating_base_;

};

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_STATE_EQUATION_HPP_ 