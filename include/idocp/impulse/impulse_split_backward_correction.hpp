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
/// @brief Backward correction for an impulse stage.
///
class ImpulseSplitBackwardCorrection {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
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

  ///
  /// @brief Coarse updates the split solution of this impulse stage. 
  /// @param[in] robot Robot model. 
  /// @param[in] aux_mat_next Auxiliary matrix of the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] s_new Coarse updated split solution of this impulse stage.
  ///
  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, 
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    ImpulseSplitKKTMatrix& kkt_matrix, 
                    const ImpulseSplitKKTResidual& kkt_residual,
                    const ImpulseSplitSolution& s, ImpulseSplitSolution& s_new);

  ///
  /// @brief Auxiliary matrix of this impulse stage. 
  /// @return const reference to the auxiliary matrix of this impulse stage. 
  ///
  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  ///
  /// @brief Performs the serial part of the backward correction. 
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in] s_new_next Coarse updated split solution of the next time stage.
  /// @param[in, out] s_new Coarse updated split solution of this impulse stage.
  ///
  void backwardCorrectionSerial(const SplitSolution& s_next, 
                                const SplitSolution& s_new_next,
                                ImpulseSplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the backward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in, out] s_new Coarse updated split solution of this impulse stage.
  ///
  void backwardCorrectionParallel(const Robot& robot, 
                                  ImpulseSplitSolution& s_new);

  ///
  /// @brief Performs the serial part of the forward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in] s_prev Split solution of the previous time stage.
  /// @param[in] s_new_prev Coarse updated split solution of the previous time 
  /// stage.
  /// @param[in, out] s_new Coarse updated split solution of this impulse stage.
  ///
  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev, 
                               const SplitSolution& s_new_prev,
                               ImpulseSplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the forward correction. 
  /// @param[in, out] s_new Coarse updated split solution of this impulse stage.
  ///
  void forwardCorrectionParallel(ImpulseSplitSolution& s_new);

  ///
  /// @brief Computes the direction from the results of the forward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] s_new Coarse updated split solution of this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  ///
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