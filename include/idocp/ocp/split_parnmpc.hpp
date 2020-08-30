#ifndef IDOCP_SPLIT_PARNMPC_HPP_ 
#define IDOCP_SPLIT_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/robot_dynamics.hpp"


namespace idocp {

///
/// @class SplitParNMPC
/// @brief Split OCP of ParNMPC for single stage. 
///
class SplitParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ///
  /// @brief Construct SplitParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr of the cost function.
  /// @param[in] cost Shared ptr of the constraints.
  ///
  SplitParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
               const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  SplitParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitParNMPC();

  ///
  /// @brief Use default copy constructor. 
  ///
  SplitParNMPC(const SplitParNMPC&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  SplitParNMPC& operator=(const SplitParNMPC&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  SplitParNMPC(SplitParNMPC&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  SplitParNMPC& operator=(SplitParNMPC&&) noexcept = default;
 
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
  /// @brief Updates the solution of the split OCP approximately. If this stage
  /// is terminal, call coarseUpdateTerminal() instead.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in] aux_mat_next_old Guess of the auxiliary matrix of the next 
  /// stage. Size must be 2*Robot::dimv() x 2*Robot::dimv().
  /// @param[out] d Split direction of this stage.
  /// @param[out] s_new_coarse Coarse updated split splition of this stage.
  ///
  void coarseUpdate(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                    const SplitSolution& s_next, 
                    const Eigen::MatrixXd& aux_mat_next_old, SplitDirection& d, 
                    SplitSolution& s_new_coarse);

  ///
  /// @brief Updates the solution of the split OCP at the terminal 
  /// approximately. If this stage is not terminal, call coarseUpdate() instead.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[out] d Split direction of this stage.
  /// @param[out] s_new_coarse Coarse updated split splition of this stage.
  ///
  void coarseUpdateTerminal(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& q_prev, 
                            const Eigen::VectorXd& v_prev, 
                            const SplitSolution& s, SplitDirection& d, 
                            SplitSolution& s_new_coarse);

  ///
  /// @brief Gets the auxiliary matrix of this stage.
  /// @param[out] aux_mat The auxiliary matrix of this stage. Size must be 
  /// 2*Robot::dimv() x 2*Robot::dimv().
  ///
  void getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const;

  ///
  /// @brief Gets the hessian of the terminal cost. Useful to initialize the 
  /// auxiliary matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[out] phixx The hessian of the terminal cost. The size must be 
  /// 2*Robot::dimv() x 2*Robot::dimv().
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
                                  SplitSolution& s_new);

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
                                 SplitSolution& s_new);

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
  void computePrimalAndDualDirection(const Robot& robot, const double dtau,
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
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage. The stage cost is recomputed. The violation of the  
  /// constriants is not computed. Instead, the previously computed residual  
  /// computed by coarseUpdate() or computeSquaredKKTErrorNorm(), is used.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolation(Robot& robot, const double t, 
                                             const double dtau, 
                                             const SplitSolution& s);

  ///
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage under step_size. The split solution of this stage and the 
  /// state of the previous stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[out] s_tmp The temporary split solution.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const Eigen::VectorXd& q_prev, 
                                             const Eigen::VectorXd& v_prev, 
                                             const SplitSolution& s, 
                                             const SplitDirection& d, 
                                             SplitSolution& s_tmp);

  ///
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage under step_size. The split solution of this stage and the 
  /// state of the previous stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s_prev Split solutio of the previous stage.
  /// @param[in] d_prev Split direction of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[out] s_tmp The temporary split solution.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const SplitSolution& s_prev, 
                                             const SplitDirection& d_prev, 
                                             const SplitSolution& s, 
                                             const SplitDirection& d, 
                                             SplitSolution& s_tmp);

  ///
  /// @brief Returns the sum of the stage cost and terminal cost, and 
  /// L1-norm of the violation of the constraints of this stage. The stage cost 
  /// and the terminal cost are recomputed. The violation of the constriants is 
  /// not computed. Instead, the previously computed residual computed by 
  /// coarseUpdateTerminal() or computeSquaredKKTErrorNormTerminal(), is used.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolationTerminal(Robot& robot, 
                                                     const double t, 
                                                     const double dtau, 
                                                     const SplitSolution& s);

  ///
  /// @brief Returns the sum of the stage cost and terminal cost and L1-norm of 
  /// the violation of constraints of this stage under step_size. The split 
  /// solution of this stage and the state of the previous stage are computed 
  /// by step_size temporary. The stage cost, the terminal cost the violation 
  /// of the constriants are computed based on the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s_prev Split solutio of the previous stage.
  /// @param[in] d_prev Split direction of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[out] s_tmp The temporary split solution.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolationTerminal(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const SplitSolution& s_prev, const SplitDirection& d_prev,
      const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp);

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] riccati Riccati factorization of this stage.
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
  /// @brief Gets the state-feedback gain for the control input torques.
  /// @param[out] Kq Gain with respec to the configuration. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Kv Gain with respec to the velocity. Size must be
  /// Robot::dimv() x Robot::dimv().
  ///
  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// residual computed by coarseUpdate(). The result is not exactly the same 
  /// as the squared KKT error norm of the original OCP. The result is the
  /// squared norm of the condensed residual. However, this variables is 
  /// sufficiently close to the original KKT error norm.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @return The squared norm of the condensed KKT residual.
  ///
  double condensedSquaredKKTErrorNorm(Robot& robot, const double t,  
                                      const double dtau, 
                                      const SplitSolution& s);

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// residual computed by coarseUpdateTerminal(). The result is not exactly  
  /// the same as the squared KKT error norm of the original OCP. The result 
  /// is the squared norm of the condensed residual. However, this variables is 
  /// sufficiently close to the original KKT error norm.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @return The squared norm of the condensed KKT residual.
  ///
  double condensedSquaredKKTErrorNormTerminal(Robot& robot, const double t, 
                                              const double dtau, 
                                              const SplitSolution& s);

  ///
  /// @brief Computes and returns the squared KKT error norm of the OCP of this
  /// stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @return The squared norm of the kKT residual.
  ///
  double computeSquaredKKTErrorNorm(Robot& robot, const double t, 
                                    const double dtau, 
                                    const Eigen::VectorXd& q_prev, 
                                    const Eigen::VectorXd& v_prev, 
                                    const SplitSolution& s,
                                    const SplitSolution& s_next);

  ///
  /// @brief Computes and returns the squared KKT error norm of the OCP of this
  /// stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage. Size must be 
  /// Robot::dimq().
  /// @param[in] v_prev Velocity of the previous stage. Size must be 
  /// Robot::dimv().
  /// @param[in] s Split solution of this stage.
  /// @return The squared norm of the kKT residual at the terminal stage.
  ///
  double computeSquaredKKTErrorNormTerminal(Robot& robot, const double t, 
                                            const double dtau, 
                                            const Eigen::VectorXd& q_prev, 
                                            const Eigen::VectorXd& v_prev, 
                                            const SplitSolution& s);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  StateEquation state_equation_;
  RobotDynamics robot_dynamics_;
  int dimv_, dimx_, dimKKT_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  Eigen::VectorXd x_res_; /// @brief Residual of state and costate used in the forward and backward correction.
  Eigen::VectorXd dx_; /// @brief Correction term of state and costate used in the forward and backward correction.
  bool use_kinematics_;

};

} // namespace idocp


#endif // IDOCP_SPLIT_PARNMPC_HPP_