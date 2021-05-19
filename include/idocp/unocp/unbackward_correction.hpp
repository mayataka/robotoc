#ifndef IDOCP_UNBACKWARD_CORRECTION_HPP_
#define IDOCP_UNBACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unbackward_correction.hpp"
#include "idocp/unocp/split_unparnmpc.hpp"
#include "idocp/unocp/terminal_unparnmpc.hpp"
#include "idocp/unocp/unconstrained_container.hpp"


namespace idocp {

///
/// @class BackwardCorrectionSolver
/// @brief Backward correction solver for optimal control problems of 
/// unconstrained rigid-body systems.
///
class UnBackwardCorrection {
public:
  ///
  /// @brief Construct a backward correction solver.
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  UnBackwardCorrection(const Robot& robot, const double T, const int N, 
                       const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  UnBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~UnBackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnBackwardCorrection(const UnBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnBackwardCorrection& operator=(const UnBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnBackwardCorrection(UnBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnBackwardCorrection& operator=(UnBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Initializes the auxiliary matrices by the terminal cost Hessian 
  /// computed by the current solution. 
  /// @param[in] robots std::vector of Robot.
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] s Solution. 
  ///
  void initAuxMat(std::vector<Robot, Eigen::aligned_allocator<Robot>>& robots, UnParNMPC& parnmpc, 
                  const double t, const UnSolution& s);

  ///
  /// @brief Linearizes the optimal control problem and coarse updates the 
  /// solution in parallel. 
  /// @param[in] robots std::vector of Robot.
  /// @param[in, out] parnmpc Optimal control problem.
  /// @param[in] t Initial time of the horizon.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in, out] unkkt_matrix KKT matrix. 
  /// @param[in, out] unkkt_residual KKT residual. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void coarseUpdate(std::vector<Robot, Eigen::aligned_allocator<Robot>>& robots, UnParNMPC& parnmpc, 
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, UnKKTMatrix& unkkt_matrix, 
                    UnKKTResidual& unkkt_residual,
                    const UnSolution& s, UnDirection& d);

  ///
  /// @brief Performs the backward correction for coarse updated solution and 
  /// computes the Newton direction. 
  /// @param[in] robots std::vector of Robot.
  /// @param[in] parnmpc Optimal control problem.
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void backwardCorrection(std::vector<Robot, Eigen::aligned_allocator<Robot>>& robots, UnParNMPC& parnmpc, 
                          const UnSolution& s, UnDirection& d);

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
  UnBackwardCorrector corrector_;
  UnSolution s_new_;
  std::vector<Eigen::MatrixXd> aux_mat_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace idocp

#endif // IDOCP_UNBACKWARD_CORRECTION_HPP_ 