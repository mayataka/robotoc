#ifndef IDOCP_TERMINAL_OCP_HPP_
#define IDOCP_TERMINAL_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

///
/// @class TerminalOCP
/// @brief Split optimal control problem at the terminal stage. 
///
class TerminalOCP {
public:
  ///
  /// @brief Construct a terminal optimal control problem.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  TerminalOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  TerminalOCP();
  
  ///
  /// @brief Destructor. 
  ///
  ~TerminalOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  TerminalOCP(const TerminalOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TerminalOCP& operator=(const TerminalOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TerminalOCP(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TerminalOCP& operator=(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  ///
  /// @brief Linearize the terminal optimal control problem around the current 
  /// solution for Newton's method, i.e, computes KKT residual and Hessian.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of the terminal stage. 
  /// @param[in] s Split solution of the terminal stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const SplitSolution& s);

  ///
  /// @brief Computes the Riccati factorization of this terminal stage. 
  /// TerminalOCP::linearizeOCP() must be called before calling this function.
  /// @param[out] riccati Riccati factorization of this terminal stage.
  /// 
  void backwardRiccatiRecursion(RiccatiSolution& riccati) const;

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this stage.
  /// @param[in] riccati Riccati factorization of this terminal stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedPrimalDirection(const RiccatiSolution& riccati,
                                       SplitDirection& d) const;

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this stage.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedDualDirection(const SplitDirection& d);

  ///
  /// @brief Returns maximum stap size of the primal variables that satisfies 
  /// the inequality constraints. TerminalOCP::computeCondensedPrimalDirection()
  /// must be called before calling this function.
  /// @return Maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  ///
  double maxPrimalStepSize();

  ///
  /// @brief Returns maximum stap size of the dual variables that satisfies 
  /// the inequality constraints. TerminalOCP::computeCondensedDualDirection()
  /// must be called before calling this function.
  /// @return Maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  ///
  double maxDualStepSize();

  ///
  /// @brief Returns the terminal cost.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  ///
  double terminalCost(Robot& robot, const double t, const SplitSolution& s);

  ///
  /// @brief Returns the terminal cost under step_size. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  ///
  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const SplitSolution& s, const SplitDirection& d);

  ///
  /// @brief Updates dual variables of the inequality constraints. 
  /// TerminalOCP::computeCondensedDualDirection() must be called before 
  /// calling this function.
  /// @param[in] step_size Dula step size of the OCP. 
  ///
  void updateDual(const double step_size);

  ///
  /// @brief Updates primal variables of this stage.
  /// TerminalOCP::computeCondensedPrimalDirection() must be called before
  /// calling this function.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double step_size, 
                    const SplitDirection& d, SplitSolution& s) const;

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  ///
  void computeKKTResidual(Robot& robot, const double t, const SplitSolution& s);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage.  
  /// SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual must be called 
  /// before calling this function.
  /// @return The squared norm of the KKT residual.
  ///
  double squaredNormKKTResidual() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  SplitSolution s_tmp_; 
  bool use_kinematics_;

};

} // namespace idocp


#endif // IDOCPT_TERMINAL_OCP_HPP_