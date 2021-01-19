#ifndef IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HPP_
#define IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hpp"
#include "idocp/impulse/impulse_split_backward_correction_data.hpp"


namespace idocp {

///
/// @class ImpulseSplitBackwardCorrection
/// @brief Split unconstrained backward correction.
///
class ImpulseSplitBackwardCorrection {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  ImpulseSplitBackwardCorrection(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitBackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitBackwardCorrection(
      const ImpulseSplitBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitBackwardCorrection& operator=(
      const ImpulseSplitBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitBackwardCorrection(
      ImpulseSplitBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitBackwardCorrection& operator=(
      ImpulseSplitBackwardCorrection&&) noexcept = default;

  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, 
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    ImpulseSplitKKTMatrix& kkt_matrix, 
                    const ImpulseSplitKKTResidual& kkt_residual,
                    const ImpulseSplitSolution& s, ImpulseSplitSolution& s_new);

  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  void backwardCorrectionSerial(const SplitSolution& s_next, 
                                const SplitSolution& s_new_next,
                                ImpulseSplitSolution& s_new);

  void backwardCorrectionParallel(const Robot& robot, 
                                  ImpulseSplitSolution& s_new);

  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev, 
                               const SplitSolution& s_new_prev,
                               ImpulseSplitSolution& s_new);

  void forwardCorrectionParallel(ImpulseSplitSolution& s_new);

  static void computeDirection(const Robot& robot, 
                               const ImpulseSplitSolution& s, 
                               const ImpulseSplitSolution& s_new, 
                               ImpulseSplitDirection& d);

private:
  int dimv_, dimx_, dimKKT_;
  ImpulseSplitKKTMatrixInverter kkt_mat_inverter_;
  ImpulseSplitBackwardCorrectionData data_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_backward_correction.hxx"

#endif // IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HPP_ 