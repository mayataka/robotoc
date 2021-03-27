#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"
#include "idocp/ocp/split_backward_correction_data.hpp"


namespace idocp {

///
/// @class SplitBackwardCorrection
/// @brief Split backward correction.
///
class SplitBackwardCorrection {
public:
  ///
  /// @brief Construct split backward correction.
  /// @param[in] robot Robot model. 
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

  ///
  /// @brief Coarse updates the split solution of this time stage. Auxiliary 
  /// matrix is assumed to zero matrix, i.e., terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time stage of time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage. 
  /// @param[in] kkt_residual Split KKT residual of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void coarseUpdate(const Robot& robot, const double dt, 
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  ///
  /// @brief Coarse updates the split solution of this time stage. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time stage of time stage. 
  /// @param[in] aux_mat_next Auxiliary matrix of the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage. 
  /// @param[in] kkt_residual Split KKT residual of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, const double dt,
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    SplitKKTMatrix& kkt_matrix, 
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
  template <typename SplitSolutionType>
  void backwardCorrectionSerial(const SplitSolutionType& s_next, 
                                const SplitSolutionType& s_new_next,
                                SplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the backward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void backwardCorrectionParallel(const Robot& robot, SplitSolution& s_new);

  ///
  /// @brief Performs the serial part of the forward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in] s_prev Split solution of the previous time stage.
  /// @param[in] s_new_prev Coarse updated split solution of the previous time 
  /// stage.
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  template <typename SplitSolutionType>
  void forwardCorrectionSerial(const Robot& robot, 
                               const SplitSolutionType& s_prev, 
                               const SplitSolutionType& s_new_prev,
                               SplitSolution& s_new);

  ///
  /// @brief Performs the parallel part of the forward correction. 
  /// @param[in, out] s_new Coarse updated split solution of this time stage.
  ///
  void forwardCorrectionParallel(SplitSolution& s_new);

  ///
  /// @brief Computes the direction from the results of the forward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_new Coarse updated split solution of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  ///
  void computeDirection(const Robot& robot, const SplitSolution& s, 
                        const SplitSolution& s_new, SplitDirection& d) const;

private:
  int dimv_, dimx_, dimu_, dimKKT_;
  bool is_switching_constraint_valid_;
  SplitKKTMatrixInverter kkt_mat_inverter_;
  SplitBackwardCorrectionData data_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp

#include "idocp/ocp/split_backward_correction.hxx"

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_ 