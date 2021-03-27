#ifndef IDOCP_BACKWARD_CORRECTION_SOLVER_HPP_
#define IDOCP_BACKWARD_CORRECTION_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"


namespace idocp {

///
/// @class BackwardCorrectionSolver
/// @brief Backward correction solver for hybrid optimal control problems.
///
class BackwardCorrectionSolver {
public:
  ///
  /// @brief Construct a backward correction solver.
  /// @param[in] robot Robot model. 
  /// @param[in] N Number of discretization of the horizon. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  BackwardCorrectionSolver(const Robot& robot, const int N, 
                           const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  BackwardCorrectionSolver();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardCorrectionSolver();
 
  ///
  /// @brief Default copy constructor. 
  ///
  BackwardCorrectionSolver(const BackwardCorrectionSolver&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  BackwardCorrectionSolver& operator=(
      const BackwardCorrectionSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardCorrectionSolver(BackwardCorrectionSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardCorrectionSolver& operator=(
      BackwardCorrectionSolver&&) noexcept = default;

  ///
  /// @brief Initializes the auxiliary matrices by the terminal cost Hessian 
  /// computed by the current solution. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  ///
  void initAuxMat(ParNMPC& parnmpc, std::vector<Robot>& robots, 
                  const Solution& s, KKTMatrix& kkt_matrix);

  ///
  /// @brief Linearizes the optimal control problem and coarse updates the 
  /// solution in parallel. 
  /// @param[in, out] parnmpc Optimal control problem.
  /// @param[in, out] corr Collection of split backward corrections.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void coarseUpdate(ParNMPC& parnmpc, BackwardCorrection& corr, 
                    std::vector<Robot>& robots, 
                    const ContactSequence& contact_sequence, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Solution& s, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual);

  ///
  /// @brief Performs the serial part of the backward correction. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in, out] corr Collection of split backward corrections.
  /// @param[in] s Solution. 
  ///
  void backwardCorrectionSerial(const ParNMPC& parnmpc, 
                                BackwardCorrection& corr, const Solution& s);

  ///
  /// @brief Performs the parallel part of the backward correction. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in, out] corr Collection of split backward corrections.
  /// @param[in] robots std::vector of Robot.
  ///
  void backwardCorrectionParallel(const ParNMPC& parnmpc, 
                                  BackwardCorrection& corr,
                                  const std::vector<Robot>& robots);

  ///
  /// @brief Performs the serial part of the forward correction. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in, out] corr Collection of split backward corrections.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] s Solution. 
  ///
  void forwardCorrectionSerial(const ParNMPC& parnmpc, BackwardCorrection& corr, 
                               const std::vector<Robot>& robots, 
                               const Solution& s);

  ///
  /// @brief Performs the parallel part of the forward correction and computes
  /// the direction. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in, out] corr Collection of split backward corrections.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void forwardCorrectionParallel(ParNMPC& parnmpc, BackwardCorrection& corr, 
                                 std::vector<Robot>& robots, 
                                 const KKTMatrix& kkt_matrix, 
                                 KKTResidual& kkt_residual, 
                                 const Solution& s, Direction& d);

  ///
  /// @brief Returns max primal step size.
  /// @return max primal step size.
  /// 
  double maxPrimalStepSize() const;

  ///
  /// @brief Returns max dual step size.
  /// @return max dual step size.
  /// 
  double maxDualStepSize() const;

  static const Eigen::VectorXd& q_prev(const ParNMPCDiscretizer& discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage < discretizer.N());
    if (discretizer.isTimeStageAfterImpulse(time_stage)) {
      return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].q;
    }
    else if (discretizer.isTimeStageAfterLift(time_stage)) {
      return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].q;
    }
    else if (time_stage > 0) {
      return s[time_stage-1].q;
    }
    else { 
      assert(time_stage == 0);
      return q;
    }
  }

  static const Eigen::VectorXd& v_prev(const ParNMPCDiscretizer& discretizer, 
                                       const Eigen::VectorXd& v, 
                                       const Solution& s, 
                                       const int time_stage) {
    assert(time_stage >= 0);
    assert(time_stage < discretizer.N());
    if (discretizer.isTimeStageAfterImpulse(time_stage)) {
      return s.impulse[discretizer.impulseIndexBeforeTimeStage(time_stage)].v;
    }
    else if (discretizer.isTimeStageAfterLift(time_stage)) {
      return s.lift[discretizer.liftIndexBeforeTimeStage(time_stage)].v;
    }
    else if (time_stage > 0) {
      return s[time_stage-1].v;
    }
    else { 
      assert(time_stage == 0);
      return v;
    }
  }

private:
  int max_num_impulse_, nthreads_, N_all_;
  Solution s_new_;
  std::vector<Eigen::MatrixXd> aux_mat_, aux_mat_impulse_, aux_mat_aux_, 
                               aux_mat_lift_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace idocp

#endif // IDOCP_BACKWARD_CORRECTION_SOLVER_HPP_ 