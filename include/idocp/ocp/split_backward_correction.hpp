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

  void coarseUpdate(const Robot& robot, const double dtau, 
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, const double dtau,
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, SplitSolution& s_new);

  template <typename MatrixType>
  void coarseUpdate(const Robot& robot, const double dtau,
                    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    SplitKKTMatrix& kkt_matrix, 
                    const SplitKKTResidual& kkt_residual,
                    const SplitSolution& s, const ImpulseSplitSolution& s_next, 
                    SplitSolution& s_new);

  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  template <typename SplitSolutionType>
  void backwardCorrectionSerial(const SplitSolutionType& s_next, 
                                const SplitSolutionType& s_new_next,
                                SplitSolution& s_new);

  void backwardCorrectionParallel(const Robot& robot, SplitSolution& s_new);

  template <typename SplitSolutionType>
  void forwardCorrectionSerial(const Robot& robot, 
                               const SplitSolutionType& s_prev, 
                               const SplitSolutionType& s_new_prev,
                               SplitSolution& s_new);

  void forwardCorrectionParallel(SplitSolution& s_new);

  static void computeDirection(const Robot& robot, const SplitSolution& s, 
                               const SplitSolution& s_new, 
                               SplitDirection& d);

  void computeDirection(const ImpulseSplitSolution& s,
                        ImpulseSplitDirection& d) const;

private:
  int dimv_, dimx_, dimu_, dimKKT_;
  bool is_impulse_condition_valid_;
  SplitKKTMatrixInverter kkt_mat_inverter_;
  SplitBackwardCorrectionData data_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp

#include "idocp/ocp/split_backward_correction.hxx"

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HPP_ 