#ifndef IDOCP_UNBACKWARD_CORRECTION_HPP_
#define IDOCP_UNBACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unbackward_correction.hpp"
#include "idocp/unocp/split_unparnmpc.hpp"
#include "idocp/unocp/terminal_unparnmpc.hpp"
#include "idocp/unocp/unconstrained_container.hpp"


namespace idocp {

///
/// @class UnBackwardCorrection
/// @brief Unconstrained backward correction.
///
class UnBackwardCorrection {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
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

  void initAuxMat(std::vector<Robot>& robots, UnParNMPC& parnmpc, 
                  const double t, const UnSolution& s);

  void coarseUpdate(std::vector<Robot>& robots, UnParNMPC& parnmpc, 
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, UnKKTMatrix& unkkt_matrix, 
                    UnKKTResidual& unkkt_residual,
                    const UnSolution& s, UnDirection& d);

  void backwardCorrection(std::vector<Robot>& robots, UnParNMPC& parnmpc, 
                          const UnSolution& s, UnDirection& d);

  double primalStepSize() const;

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