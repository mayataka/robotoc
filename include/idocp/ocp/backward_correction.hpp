#ifndef IDOCP_BACKWARD_CORRECTION_HPP_
#define IDOCP_BACKWARD_CORRECTION_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"

#include "Eigen/Core"

namespace idocp {

class BackwardCorrection {
public:
  BackwardCorrection(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  BackwardCorrection();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardCorrection();

  ///
  /// @brief Default copy constructor. 
  ///
  BackwardCorrection(const BackwardCorrection&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  BackwardCorrection& operator=(const BackwardCorrection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardCorrection(BackwardCorrection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardCorrection& operator=(BackwardCorrection&&) noexcept = default;

  void coarseUpdate(const Robot& robot, const SplitSolution& s, 
                    SplitDirection& d, KKTMatrix& kkt_matrix, 
                    const KKTResidual& kkt_residual,
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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd x_res_, dx_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  int dimv_, dimx_, dimKKT_;

};

} // namespace idocp

#include "idocp/ocp/backward_correction.hxx"

#endif // IDOCP_BACKWARD_CORRECTION_HPP_ 