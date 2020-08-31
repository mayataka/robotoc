#ifndef IDOCP_TERMINAL_OCP_HPP_
#define IDOCP_TERMINAL_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

///
/// @class TerminalOCP
/// @brief Split OCP for terminal stage. 
///
class TerminalOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ///
  /// @brief Construct a terminal OCP.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr of the cost function.
  /// @param[in] cost Shared ptr of the constraints.
  ///
  TerminalOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  TerminalOCP();
  
  ///
  /// @brief Destructor. 
  ///
  ~TerminalOCP();

  ///
  /// @brief Use default copy constructor. 
  ///
  TerminalOCP(const TerminalOCP&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  TerminalOCP& operator=(const TerminalOCP&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  TerminalOCP(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  TerminalOCP& operator=(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(const Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  ///
  /// @brief Linearize the terminal OCP for Newton's method around the current 
  /// solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[out] riccati Riccati factorization of this stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const SplitSolution& s, 
                    RiccatiFactorization& riccati);

  ///
  /// @brief Computes the Newton direction of the condensed variables of this 
  /// stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 SplitDirection& d);

  ///
  /// @brief Returns maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  ///
  double maxPrimalStepSize();

  ///
  /// @brief Returns maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  ///
  double maxDualStepSize();

  ///
  /// @brief Returns the terminal cost.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  ///
  double terminalCost(Robot& robot, const double t, const SplitSolution& s);

  ///
  /// @brief Returns the terminal cost under step_size. The split solution of 
  /// this stage and is computed by step_size temporary. 
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  ///
  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const SplitSolution& s, const SplitDirection& d);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] step_size Dula step size of the OCP. 
  ///
  void updateDual(const double step_size);

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] riccati Riccati factorization of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double step_size, 
                    const RiccatiFactorization& riccati,
                    const SplitDirection& d, SplitSolution& s) const;

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// KKT residual computed by linearizeOCP().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @return The squared norm of the KKT residual.
  ///
  double squaredKKTErrorNorm(Robot& robot, const double t, 
                             const SplitSolution& s) const;

  ///
  /// @brief Computes and returns the squared KKT error norm by using  
  /// previously computed KKT residual computed by linearizeOCP().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @return The squared norm of the KKT residual.
  ///
  double computeSquaredKKTErrorNorm(Robot& robot, const double t, 
                                    const SplitSolution& s);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  SplitSolution s_tmp_; /// @brief Temporary split solution used in line search.
  bool use_kinematics_;

};

} // namespace idocp


#endif // IDOCPT_TERMINAL_OCP_HPP_