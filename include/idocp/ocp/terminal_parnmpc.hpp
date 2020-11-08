#ifndef IDOCP_TERMINAL_PARNMPC_HPP_ 
#define IDOCP_TERMINAL_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/state_equation.hpp"
// #include "idocp/ocp/robot_dynamics.hpp"
#include "idocp/ocp/contact_dynamics.hpp"
#include "idocp/ocp/backward_correction.hpp"


namespace idocp {

///
/// @class TerminalParNMPC
/// @brief Split optimal control problem of ParNMPC algorithm for the terminal 
/// stage. 
///
class TerminalParNMPC {
public:
  ///
  /// @brief Construct TerminalParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  TerminalParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                  const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  TerminalParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~TerminalParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  TerminalParNMPC(const TerminalParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TerminalParNMPC& operator=(const TerminalParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TerminalParNMPC(TerminalParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TerminalParNMPC& operator=(TerminalParNMPC&&) noexcept = default;

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
  void initConstraints(Robot& robot, const int time_step, const double dtau, 
                       const SplitSolution& s);

  ///
  /// @brief Linearize the terminal OCP. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  ///
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status,
                    const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, const SplitSolution& s);

  ///
  /// @brief Updates the solution of the split OCP approximately. Call  
  /// linearizeOCP() or linearizeTerminalOCP() before calling this function.
  /// @param[in] aux_mat_next_old Guess of the auxiliary matrix of the next 
  /// stage. Size must be 2 * Robot::dimv() x 2 * Robot::dimv().
  /// @param[out] d Split direction of this stage.
  /// @param[out] s_new_coarse Coarse updated split splition of this stage.
  ///
  void coarseUpdate(SplitDirection& d, SplitSolution& s_new_coarse);

  ///
  /// @brief Gets the auxiliary matrix of this stage.
  /// @param[out] aux_mat The auxiliary matrix of this stage. Size must be 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  void getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const;

  ///
  /// @brief Gets the hessian of the terminal cost. Useful to initialize the 
  /// auxiliary matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[out] phixx The hessian of the terminal cost. The size must be 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  void getTerminalCostHessian(Robot& robot, const double t, 
                              const SplitSolution& s, Eigen::MatrixXd& phixx);

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
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  ///
  /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  /// after backwardCorrectionSerial() and before calling 
  /// forwardCorrectionSerial().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] d Split directtion of the current stage.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void backwardCorrectionParallel(const Robot& robot, SplitDirection& d,
                                  SplitSolution& s_new) const;

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
                               const SplitSolution& s_new_prev, 
                               SplitSolution& s_new);

  ///
  /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  /// after forwardCorrectionSerial() and before calling 
  /// computePrimalAndDualDirection().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] d Split directtion of the current stage.
  /// @param[out] s_new Split solution of the current stage at the current 
  /// iteration.
  ///
  void forwardCorrectionParallel(const Robot& robot, SplitDirection& d, 
                                 SplitSolution& s_new) const;

  ///
  /// @brief Computes the direction of the primal and dual solution. Call after 
  /// forwardCorrectionParallel() and before calling 
  /// computePrimalAndDualDirection().
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dtau Discretization length of the OCP. Must be positive.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_new Corrected split solution of the current stage.
  /// @param[out] d Split direction of this stage.
  /// 
  void computePrimalAndDualDirection(Robot& robot, const double dtau,
                                     const SplitSolution& s,
                                     const SplitSolution& s_new,
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
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const SplitDirection& d, SplitSolution& s);

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
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  ///
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                          const double t, const double dtau, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, 
                          const SplitSolution& s);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const double dtau) const;

  ///
  /// @brief Computes the stage cost of this stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] primal_step_size Primal step size of the OCP. Default is 0.
  /// @return Stage cost of this stage.
  /// 
  double terminalCost(Robot& robot, const double t, const double dtau, 
                      const SplitSolution& s, const double primal_step_size=0);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_next Configuration at the next stage.
  /// @param[in] v_next Generaized velocity at the next stage.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dtau, 
                             const SplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

  ///
  /// @brief Gets the state-feedback gain for the control input torques.
  /// @param[out] Kq Gain with respec to the configuration. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Kv Gain with respec to the velocity. Size must be
  /// Robot::dimv() x Robot::dimv().
  ///
  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  // RobotDynamics robot_dynamics_;
  ContactDynamics contact_dynamics_;
  BackwardCorrection backward_correction_;
  int dimv_, dimx_, dimKKT_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  Eigen::VectorXd x_res_; /// @brief Residual of state and costate used in the forward and backward correction.
  Eigen::VectorXd dx_; /// @brief Correction term of state and costate used in the forward and backward correction.
  bool use_kinematics_;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  inline void setContactStatusForKKT(const ContactStatus& contact_status) {
    kkt_residual_.setContactStatus(contact_status);
    kkt_matrix_.setContactStatus(contact_status);
  }

};

} // namespace idocp

#endif // IDOCP_TERMINAL_PARNMPC_HPP_ 