#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"

#include "Eigen/Core"

namespace idocp {

class SplitBackwardCorrection {
public:
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
  /// @brief Default copy assign operator. 
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

  void coarseUpdate(const Robot& robot, const SplitSolution& s, 
                    SplitDirection& d, SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    SplitSolution& s_new_coarse);

  template <typename SplitSolutionType>
  void backwardCorrectionSerial(const SplitSolutionType& s_next,
                                const SplitSolutionType& s_new_next,
                                SplitSolution& s_new);

  void backwardCorrectionParallel(const Robot& robot, SplitDirection& d,
                                  SplitSolution& s_new) const;

  template <typename SplitSolutionType>
  void forwardCorrectionSerial(const Robot& robot, 
                               const SplitSolutionType& s_prev,
                               const SplitSolutionType& s_new_prev, 
                               SplitSolution& s_new);

  void forwardCorrectionParallel(SplitDirection& d, 
                                 SplitSolution& s_new) const;

  static void computeDirection(const Robot& robot, const SplitSolution& s, 
                               const SplitSolution& s_new, SplitDirection& d);

  template <typename MatrixType>
  void getAuxiliaryMatrix(const Eigen::MatrixBase<MatrixType>& aux_mat) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd x_res_, dx_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  int dimv_, dimx_, dimKKT_;

};

} // namespace idocp

#include "idocp/ocp/split_backward_correction.hxx"

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_ 