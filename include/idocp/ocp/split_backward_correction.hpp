#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"


namespace idocp {

///
/// @class SplitBackwardCorrection
/// @brief Split unconstrained backward correction.
///
class SplitBackwardCorrection {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  SplitBackwardCorrection(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~SplitBackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitBackwardCorrection(const SplitBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitBackwardCorrection& operator=(const SplitBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitBackwardCorrection(SplitBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitBackwardCorrection& operator=(SplitBackwardCorrection&&) noexcept = default;

  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, const double dtau,
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitDirection& d, 
                    SplitSolution& s_new);

  void coarseUpdate(const Robot& robot, const double dtau, 
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitDirection& d, 
                    SplitSolution& s_new);

  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  void backwardCorrectionSerial(const SplitSolution& s_next, 
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  void backwardCorrectionParallel(const Robot& robot, SplitDirection& d, 
                                  SplitSolution& s_new) const;

  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev, 
                               const SplitSolution& s_new_prev,
                               SplitSolution& s_new);

  void forwardCorrectionParallel(SplitDirection& d, SplitSolution& s_new) const;

  static void computeDirection(const Robot& robot, const SplitSolution& s, 
                               const SplitSolution& s_new, SplitDirection& d);

private:
  int dimv_, dimx_, dimKKT_;
  SplitKKTMatrixInverter kkt_mat_inverter_;
  Eigen::MatrixXd KKT_mat_inv_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp

#include "idocp/ocp/split_backward_correction.hxx"

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_ 