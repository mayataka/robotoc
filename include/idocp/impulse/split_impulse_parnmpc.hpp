#ifndef IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_ 
#define IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"


namespace idocp {

///
/// @class SplitParNMPC
/// @brief Split optimal control problem of ParNMPC algorithm for single stage. 
///
class SplitImpulseParNMPC {
public:
  ///
  /// @brief Construct SplitParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitImpulseParNMPC(const Robot& robot, 
                      const std::shared_ptr<ImpulseCostFunction>& cost,
                      const std::shared_ptr<ImpulseConstraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  SplitImpulseParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitImpulseParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitImpulseParNMPC(const SplitImpulseParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitImpulseParNMPC& operator=(const SplitImpulseParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitImpulseParNMPC(SplitImpulseParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitImpulseParNMPC& operator=(SplitImpulseParNMPC&&) noexcept = default;
 
  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Updates the solution of the split OCP approximately. If this stage
  /// is terminal, call coarseUpdateTerminal() instead.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in] aux_mat_next_old Guess of the auxiliary matrix of the next 
  /// stage. Size must be 2 * Robot::dimv() x 2 * Robot::dimv().
  /// @param[out] d Split direction of this stage.
  /// @param[out] s_new_coarse Coarse updated split splition of this stage.
  ///
  void coarseUpdate(Robot& robot, const ContactStatus& contact_status,
                    const double t, const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, 
                    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
                    const Eigen::MatrixXd& aux_mat_next_old, 
                    ImpulseSplitDirection& d, 
                    ImpulseSplitSolution& s_new_coarse);

  ///
  /// @brief Gets the auxiliary matrix of this stage.
  /// @param[out] aux_mat The auxiliary matrix of this stage. Size must be 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  void getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const;

  ///
  /// @brief Corrects the part of the solution updated coarsely. Call serially 
  /// after calling coarseUpdate() or coarseUpdateTerminal() and before calling
  /// backwardCorrectionParallel().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s_next Split solution of the next stage at the previous 
  /// iteration.
  /// @param[in] s_new_next Split solution of the next stage at the current 
  /// iteration.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void backwardCorrectionSerial(const Robot& robot, const SplitSolution& s_next,
                                const SplitSolution& s_next_new,
                                ImpulseSplitSolution& s_new);

  ///
  /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  /// after backwardCorrectionSerial() and before calling 
  /// forwardCorrectionSerial().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] d Split directtion of the current stage.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void backwardCorrectionParallel(const Robot& robot, ImpulseSplitDirection& d,
                                  ImpulseSplitSolution& s_new) const;

  ///
  /// @brief Corrects the part of the solutio updated coarsely. Call serially 
  /// after backwardCorrectionParallel() and before calling 
  /// forwardCorrectionParallel().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s_prev Split solution of the previous stage at the previous 
  /// iteration.
  /// @param[in] s_new_prev Split solution of the previous stage at the current 
  /// iteration.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev,
                               const SplitSolution& s_prev_new, 
                               ImpulseSplitSolution& s_new);

  ///
  /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  /// after forwardCorrectionSerial() and before calling 
  /// computePrimalAndDualDirection().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] d Split directtion of the current stage.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void forwardCorrectionParallel(const Robot& robot, ImpulseSplitDirection& d, 
                                 ImpulseSplitSolution& s_new) const;

  ///
  /// @brief Computes the direction of the primal and dual solution. Call after 
  /// forwardCorrectionParallel() and before calling 
  /// computePrimalAndDualDirection().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_new Corrected split solution of the current stage.
  /// @param[out] d Split direction of this stage.
  /// 
  void computePrimalAndDualDirection(Robot& robot, 
                                     const ImpulseSplitSolution& s,
                                     const ImpulseSplitSolution& s_new,
                                     ImpulseSplitDirection& d);

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
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage. The stage cost is recomputed. The violation of the  
  /// constriants is not computed. Instead, the previously computed residual  
  /// computed by coarseUpdate() or computeSquaredKKTErrorNorm(), is used.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndConstraintViolation(
      Robot& robot, const double t, const ImpulseSplitSolution& s);

  ///
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage under step_size. The split solution of this stage and the 
  /// state of the previous stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[out] s_tmp The temporary split solution.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndConstraintViolation(
      Robot& robot, const ContactStatus& contact_status, const double step_size, 
      const double t, const Eigen::VectorXd& q_prev, 
      const Eigen::VectorXd& v_prev, const ImpulseSplitSolution& s, 
      const ImpulseSplitDirection& d, ImpulseSplitSolution& s_tmp);

  ///
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage under step_size. The split solution of this stage and the 
  /// state of the previous stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s_prev Split solutio of the previous stage.
  /// @param[in] d_prev Split direction of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[out] s_tmp The temporary split solution.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const SplitSolution& s_prev, const SplitDirection& d_prev, 
    const ImpulseSplitSolution& s, const ImpulseSplitDirection& d, 
    ImpulseSplitSolution& s_tmp);

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double step_size, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] step_size Dula step size of the OCP. 
  ///
  void updateDual(const double step_size);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status,
                          const double t, const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, 
                          const ImpulseSplitSolution& s, 
                          const SplitSolution& s_next);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<ImpulseCostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<ImpulseConstraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseKKTResidual kkt_residual_;
  ImpulseKKTMatrix kkt_matrix_;
  ImpulseDynamicsBackwardEuler impulse_dynamics_;
  int dimv_, dimx_, dimKKT_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  Eigen::VectorXd x_res_; /// @brief Residual of state and costate used in the forward and backward correction.
  Eigen::VectorXd dx_; /// @brief Correction term of state and costate used in the forward and backward correction.

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  inline void setContactStatusForKKT(const ContactStatus& contact_status) {
    kkt_residual_.setContactStatus(contact_status);
    kkt_matrix_.setContactStatus(contact_status);
  }

  double cost(Robot& robot, const double t, const ImpulseSplitSolution& s);

  double constraintViolation() const;

};

} // namespace idocp


#endif // IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_ 