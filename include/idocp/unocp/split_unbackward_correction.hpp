#ifndef IDOCP_SPLIT_UNBACKWARD_CORRECTION_HPP_
#define IDOCP_SPLIT_UNBACKWARD_CORRECTION_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unkkt_matrix_inverter.hpp"


namespace idocp {

///
/// @class SplitUnBackwardCorrection
/// @brief Split unconstrained backward correction.
///
class SplitUnBackwardCorrection {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  SplitUnBackwardCorrection(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitUnBackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnBackwardCorrection();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnBackwardCorrection(const SplitUnBackwardCorrection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitUnBackwardCorrection& operator=(const SplitUnBackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitUnBackwardCorrection(SplitUnBackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnBackwardCorrection& operator=(SplitUnBackwardCorrection&&) noexcept = default;

  template <typename MatrixType>
  void coarseUpdate(const Eigen::MatrixBase<MatrixType>& aux_mat_next,
                    const double dt, SplitUnKKTMatrix& unkkt_matrix, 
                    const SplitUnKKTResidual& unkkt_residual,
                    const SplitSolution& s, SplitDirection& d, 
                    SplitSolution& s_new);

  void coarseUpdate(const double dt, SplitUnKKTMatrix& unkkt_matrix, 
                    const SplitUnKKTResidual& unkkt_residual,
                    const SplitSolution& s, SplitDirection& d, 
                    SplitSolution& s_new);

  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  void backwardCorrectionSerial(const SplitSolution& s_next, 
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  void backwardCorrectionParallel(SplitDirection& d, 
                                  SplitSolution& s_new) const;

  void forwardCorrectionSerial(const SplitSolution& s_prev, 
                               const SplitSolution& s_new_prev,
                               SplitSolution& s_new);

  void forwardCorrectionParallel(SplitDirection& d, SplitSolution& s_new) const;

  static void computeDirection(const SplitSolution& s, 
                               const SplitSolution& s_new, 
                               SplitDirection& d);

private:
  int dimv_, dimx_, dimKKT_;
  SplitUnKKTMatrixInverter kkt_mat_inverter_;
  Eigen::MatrixXd KKT_mat_inv_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp

#include "idocp/unocp/split_unbackward_correction.hxx"

#endif // IDOCP_SPLIT_UNBACKWARD_CORRECTION_HPP_ 