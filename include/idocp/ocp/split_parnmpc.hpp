#ifndef IDOCP_SPLIT_PARNMPC_HPP_ 
#define IDOCP_SPLIT_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

///
/// @class SplitParNMPC
/// @brief Split optimal control problem of ParNMPC algorithm for single stage. 
///
class SplitParNMPC {
public:
  ///
  /// @brief Construct SplitParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
               const std::shared_ptr<Constraints>& constraints, 
               const double dtau);

  ///
  /// @brief Default constructor. 
  ///
  SplitParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitParNMPC(const SplitParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitParNMPC& operator=(const SplitParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitParNMPC(SplitParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitParNMPC& operator=(SplitParNMPC&&) noexcept = default;

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
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const int time_step, 
                       const SplitSolution& s);


  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  ///
  template <typename SplitSolutionType>
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, 
                    const SplitSolution& s, const SplitSolutionType& s_next, 
                    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian. Also linearize the 
  /// switching constraints.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status at this stage. 
  /// @param[in] impulse_status Impulse status of the next impulse. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  ///
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const ImpulseStatus& impulse_status, const double t, 
                    const double dtau, const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                    const ImpulseSplitSolution& s_next, 
                    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables  
  /// at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot, const double dtau, 
                                       const SplitSolution& s, 
                                       SplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables 
  /// at this stage.
  /// @tparam SplitDirectionType Type of the split direction at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] kkt_matrix KKT matrix of this stage.
  /// @param[in] kkt_residual KKT residual of this stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedDualDirection(const Robot& robot, const double dtau, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual,
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
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] primal_step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size of the OCP. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  template <typename SplitSolutionType>
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                          const double t, const double dtau, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                          const SplitSolutionType& s_next, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status at this stage. 
  /// @param[in] impulse_status Impulse status of the next impulse. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                          const ImpulseStatus& impulse_status, const double t, 
                          const double dtau, const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                          const ImpulseSplitSolution& s_next, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const SplitKKTResidual& kkt_residual, 
                                const double dtau) const;

  ///
  /// @brief Computes the stage cost of this stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] primal_step_size Primal step size of the OCP. Default is 0.
  /// @return Stage cost of this stage.
  /// 
  double stageCost(Robot& robot, const double t, const double dtau, 
                   const SplitSolution& s, const double primal_step_size=0);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] kkt_residual KKT residual of this stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev,
                             const Eigen::VectorXd& v_prev,
                             const SplitSolution& s,
                             SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] impulse_status Impulse status of the next impulse. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] kkt_residual KKT residual of this stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const ImpulseStatus& impulse_status, 
                             const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev,
                             const Eigen::VectorXd& v_prev,
                             const SplitSolution& s,
                             SplitKKTResidual& kkt_residual);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ContactDynamics contact_dynamics_;
  bool use_kinematics_, has_floating_base_;
  double stage_cost_, constraint_violation_;

};

} // namespace idocp

#include "idocp/ocp/split_parnmpc.hxx"

#endif // IDOCP_SPLIT_PARNMPC_HPP_