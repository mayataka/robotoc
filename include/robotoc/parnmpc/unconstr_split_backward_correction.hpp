#ifndef ROBOTOC_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_
#define ROBOTOC_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/parnmpc/unconstr_kkt_matrix_inverter.hpp"


namespace robotoc {

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
  void coarseUpdate(const Eigen::MatrixXd& aux_mat_next,
                    const double dt, const SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  ///
  /// @brief Coarse updates the split solution of this time stage. 
  /// @param[in] aux_mat_next Auxiliary matrix of the next time stage. 
  /// @param[in] dt Time stage of time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage. 
  /// @param[in] kkt_residual Split KKT residual of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void coarseUpdate(const Eigen::Block<Eigen::MatrixXd>& aux_mat_next,
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


  template <typename MatrixType>
  void coarseUpdate_impl(const Eigen::MatrixBase<MatrixType>& aux_mat_next, 
                         const double dt, const SplitKKTMatrix& kkt_matrix, 
                         const SplitKKTResidual& kkt_residual, 
                         const SplitSolution& s, SplitSolution& s_new) {
    assert(aux_mat_next.rows() == dimx_);
    assert(aux_mat_next.cols() == dimx_);
    H_.topLeftCorner(dimv_, dimv_)     = kkt_matrix.Qaa;
    H_.bottomLeftCorner(dimx_, dimv_)  = kkt_matrix.Qxu; // This is actually Qxa
    H_.bottomRightCorner(dimx_, dimx_) = kkt_matrix.Qxx.transpose();
    H_.bottomRightCorner(dimx_, dimx_).noalias() += aux_mat_next;
    kkt_mat_inverter_.invert(dt, H_, kkt_mat_inv_);
    kkt_res_.head(dimx_)           = kkt_residual.Fx;
    kkt_res_.segment(dimx_, dimv_) = kkt_residual.la;
    kkt_res_.tail(dimx_)           = kkt_residual.lx;
    d_.noalias() = kkt_mat_inv_ * kkt_res_;
    s_new.lmd = s.lmd - d_.head(dimv_);
    s_new.gmm = s.gmm - d_.segment(dimv_, dimv_);
    s_new.a   = s.a   - d_.segment(2*dimv_, dimv_); 
    s_new.q   = s.q   - d_.segment(3*dimv_, dimv_);
    s_new.v   = s.v   - d_.tail(dimv_);
  }

};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HPP_ 