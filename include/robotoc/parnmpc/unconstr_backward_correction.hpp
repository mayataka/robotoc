#ifndef ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_
#define ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/parnmpc/unconstr_split_backward_correction.hpp"
#include "robotoc/unconstr/unconstr_ocp_data.hpp"
#include "robotoc/unconstr/parnmpc_intermediate_stage.hpp"
#include "robotoc/unconstr/parnmpc_terminal_stage.hpp"


namespace robotoc {

///
/// @class UnconstrBackwardCorrection
/// @brief Backward correction for optimal control problems of 
/// unconstrained rigid-body systems.
///
class UnconstrBackwardCorrection {
public:
  ///
  /// @brief Construct a backward correction.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  UnconstrBackwardCorrection(const OCP& ocp, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrBackwardCorrection() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrBackwardCorrection(const UnconstrBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrBackwardCorrection& operator=(const UnconstrBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrBackwardCorrection(UnconstrBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrBackwardCorrection& operator=(UnconstrBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Initializes the auxiliary matrices by the terminal cost Hessian 
  /// computed by the current solution. 
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void initAuxMat(aligned_vector<Robot>& robots, 
                  const std::vector<GridInfo>& time_discretization,
                  const Solution& s, KKTMatrix& kkt_matrix, 
                  KKTResidual& kkt_residual);

  ///
  /// @brief Computes the cost and constraint violations. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(aligned_vector<Robot>& robots, 
                       const std::vector<GridInfo>& time_discretization, 
                       const Solution& s);

  ///
  /// @brief Computes the cost and constraint violations. 
  /// @param[in, out] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalOCP(aligned_vector<Robot>& robots, 
               const std::vector<GridInfo>& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual and matrix. 
  /// @param[in, out] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalKKT(aligned_vector<Robot>& robots, 
               const std::vector<GridInfo>& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTMatrix& kkt_matrix, 
               KKTResidual& kkt_residual);

  ///
  /// @brief Eval KKT and coarse updates the solution leveraging the parallel
  /// computation. 
  /// @param[in, out] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void coarseUpdate(aligned_vector<Robot>& robots, 
                    const std::vector<GridInfo>& time_discretization,
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Solution& s, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual);

  ///
  /// @brief Performs the backward correction for coarse updated solution and 
  /// computes the Newton direction. 
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void backwardCorrection(const std::vector<GridInfo>& time_discretization,
                          const Solution& s, const KKTMatrix& kkt_matrix, 
                          const KKTResidual& kkt_residual, Direction& d);

  ///
  /// @brief Gets the maximum primal step size of the fraction-to-boundary-rule.
  /// @return The primal step size of the fraction-to-boundary-rule.
  ///
  double primalStepSize() const;

  ///
  /// @brief Gets the maximum dual step size of the fraction-to-boundary-rule.
  /// @return The dual step size of the fraction-to-boundary-rule.
  ///
  double dualStepSize() const;

  ///
  /// @brief Gets the performance index of the evaluation. 
  /// @return const reference to the performance index.
  ///
  const PerformanceIndex& getEval() const {
    return performance_index_;
  }

  ///
  /// @brief Integrates the solution. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(const aligned_vector<Robot>& robots,
                         const std::vector<GridInfo>& time_discretization, 
                         const double primal_step_size, 
                         const double dual_step_size, Direction& d, Solution& s);

private:
  int nthreads_;
  OCP ocp_;
  ParNMPCIntermediateStage intermediate_stage_;
  ParNMPCTerminalStage terminal_stage_;
  aligned_vector<UnconstrOCPData> data_;
  PerformanceIndex performance_index_;
  aligned_vector<UnconstrSplitBackwardCorrection> corrector_;
  Solution s_new_;
  std::vector<Eigen::MatrixXd> aux_mat_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_ 