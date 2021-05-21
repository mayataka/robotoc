#ifndef IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_
#define IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/parnmpc/unconstr_kkt_matrix_inverter.hpp"


namespace idocp {

///
/// @class UnconstrSplitBackwardCorrection 
/// @brief Split backward correction.
///
class UnconstrSplitBackwardCorrection {
public:
  ///
  /// @brief Construct split backward correction.
  /// @param[in] robot Robot model. 
  ///
  UnconstrSplitBackwardCorrection(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrSplitBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrSplitBackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrSplitBackwardCorrection(
      const UnconstrSplitBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrSplitBackwardCorrection& operator=(
      const UnconstrSplitBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrSplitBackwardCorrection(
      UnconstrSplitBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrSplitBackwardCorrection& operator=(
      UnconstrSplitBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Coarse updates the split solution of this time stage. 
  /// @param[in] aux_mat_next Auxiliary matrix of the next time stage. 
  /// @param[in] dt Time stage of time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage. 
  /// @param[in] kkt_residual Split KKT residual of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  template <typename MatrixType>
  void coarseUpdate(const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    const double dt, const SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  ///
  /// @brief Coarse updates the split solution of this time stage. Auxiliary 
  /// matrix is assumed to zero matrix, i.e., terminal stage.
  /// @param[in] dt Time stage of time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage. 
  /// @param[in] kkt_residual Split KKT residual of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void coarseUpdate(const double dt, const SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  ///
  /// @brief Auxiliary matrix of this time stage. 
  /// @return const reference to the auxiliary matrix of this time stage. 
  ///
  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  ///
  /// @brief Performs the serial part of the backward correction. 
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in] s_new_next Coarse updated split solution of the next time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void backwardCorrectionSerial(const SplitSolution& s_next, 
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the backward correction. 
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void backwardCorrectionParallel(SplitSolution& s_new);

  ///
  /// @brief Performs the serial part of the forward correction. 
  /// @param[in] s_prev Split solution of the previous time stage.
  /// @param[in] s_new_prev Coarse updated split solution of the previous time 
  /// stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void forwardCorrectionSerial(const SplitSolution& s_prev, 
                               const SplitSolution& s_new_prev,
                               SplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the forward correction. 
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void forwardCorrectionParallel(SplitSolution& s_new);

  ///
  /// @brief Computes the direction from the results of the forward correction. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_new Coarse updated split solution of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  ///
  static void computeDirection(const SplitSolution& s, 
                               const SplitSolution& s_new, 
                               SplitDirection& d);

private:
  int dimv_, dimx_, dimkkt_;
  UnconstrKKTMatrixInverter kkt_mat_inverter_;
  Eigen::MatrixXd H_, kkt_mat_inv_;
  Eigen::VectorXd kkt_res_, d_, x_res_, dx_;

};

} // namespace idocp

#include "idocp/parnmpc/unconstr_split_backward_correction.hxx"

#endif // IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_ 