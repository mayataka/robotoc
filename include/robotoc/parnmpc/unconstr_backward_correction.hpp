#ifndef ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_
#define ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/solution.hpp"
#include "robotoc/ocp/direction.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/unconstr/unconstr_parnmpc.hpp"
#include "robotoc/parnmpc/unconstr_split_backward_correction.hpp"


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
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  UnconstrBackwardCorrection(const Robot& robot, const double T, const int N, 
                             const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrBackwardCorrection();
 
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
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void initAuxMat(aligned_vector<Robot>& robots, UnconstrParNMPC& parnmpc, 
                  const double t, const Solution& s, 
                  KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the optimal control problem and coarse updates the 
  /// solution in parallel. 
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in, out] parnmpc Optimal control problem.
  /// @param[in] t Initial time of the horizon.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in] s Solution. 
  ///
  void coarseUpdate(aligned_vector<Robot>& robots, UnconstrParNMPC& parnmpc, 
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual, const Solution& s);

  ///
  /// @brief Performs the backward correction for coarse updated solution and 
  /// computes the Newton direction. 
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in] s Solution. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void backwardCorrection(UnconstrParNMPC& parnmpc, const Solution& s, 
                          const KKTMatrix& kkt_matrix, 
                          const KKTResidual& kkt_residual, Direction& d);

  ///
  /// @brief Returns max primal step size.
  /// @return max primal step size.
  /// 
  double primalStepSize() const;

  ///
  /// @brief Returns max dual step size.
  /// @return max dual step size.
  /// 
  double dualStepSize() const;

private:
  int N_, nthreads_;
  double T_, dt_;
  std::vector<UnconstrSplitBackwardCorrection> corrector_;
  Solution s_new_;
  std::vector<Eigen::MatrixXd> aux_mat_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_BACKWARD_CORRECTION_HPP_ 