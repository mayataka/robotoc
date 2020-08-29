#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

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
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

///
/// @class SplitOCP
/// @brief Split OCP for single stage. 
///
class SplitOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ///
  /// @brief Construct a split OCP.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr of the cost function.
  /// @param[in] cost Shared ptr of the constraints.
  ///
  SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
           const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  SplitOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitOCP();

  ///
  /// @brief Use default copy constructor. 
  ///
  SplitOCP(const SplitOCP&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  SplitOCP& operator=(const SplitOCP&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  SplitOCP(SplitOCP&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  SplitOCP& operator=(SplitOCP&&) noexcept = default;

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
  /// @brief Linearize the OCP for Newton's method around the current solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                    const SplitSolution& s_next);

  ///
  /// @brief Computes the Riccati factorization of this stage from the 
  /// factorization of the previous stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] riccati_next Riccati factorization of the next stage.
  /// @param[out] riccati Riccati factorization of this stage.
  /// 
  void backwardRiccatiRecursion(const double dtau, 
                                const RiccatiFactorization& riccati_next,
                                RiccatiFactorization& riccati);

  ///
  /// @brief Computes the Newton direction of the state of this stage from the 
  /// one of the previous stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] d Split direction of this stage.
  /// @param[in] d_next Split direction of the next stage.
  /// 
  void forwardRiccatiRecursion(const double dtau, SplitDirection& d,   
                               SplitDirection& d_next);

  ///
  /// @brief Computes the Newton direction of the condensed variables of this 
  /// stage.
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedDirection(const Robot& robot, const double dtau, 
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
  /// computed by linearizeOCP() or computeSquaredKKTErrorNorm(), is used.
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
  /// state of the next stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in] d_next Split direction of the next stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const SplitSolution& s, 
                                             const SplitDirection& d,
                                             const SplitSolution& s_next, 
                                             const SplitDirection& d_next);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] step_size Dula step size of the OCP. 
  ///
  void updateDual(const double step_size);

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
                    const RiccatiFactorization& riccati, 
                    const SplitDirection& d, SplitSolution& s);

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
  /// KKT residual computed by linearizeOCP(). The result is not exactly the 
  /// same as the squared KKT error norm of the original OCP. The result is the
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
  /// @brief Computes and returns the squared KKT error norm of the OCP of this
  /// stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @return The squared norm of the kKT residual.
  ///
  double computeSquaredKKTErrorNorm(Robot& robot, const double t, 
                                    const double dtau, 
                                    const Eigen::VectorXd& q_prev, 
                                    const SplitSolution& s, 
                                    const SplitSolution& s_next);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  StateEquation state_equation_;
  RobotDynamics robot_dynamics_;
  RiccatiGain riccati_gain_;
  RiccatiMatrixFactorizer riccati_factorizer_;
  RiccatiMatrixInverter riccati_inverter_;
  Eigen::MatrixXd Ginv_; /// @brief Inverse of the Riccati matrix G.
  SplitSolution s_tmp_; /// @brief Temporary split solution used in line search.
  int dimv_, dimf_, dimc_;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = robot.dim_passive() + robot.dimf();
  }

  ///
  /// @brief Gets the block matrix of Ginv wth appropriate size. Before calling 
  /// this function, call setContactStatus() to update contact dimensions.
  /// @return Block matrix of Ginv wth appropriate size.
  ///
  inline Eigen::Block<Eigen::MatrixXd> Ginv_active() {
    return Ginv_.topLeftCorner(dimv_+dimf_+dimc_, dimv_+dimf_+dimc_);
  }

};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_