#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HPP_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitBackwardCorrectionData
/// @brief Split unconstrained backward correction.
///
class SplitBackwardCorrectionData {
public:
  ///
  /// @brief Constructor.
  ///
  SplitBackwardCorrectionData(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitBackwardCorrectionData();

  ///
  /// @brief Destructor. 
  ///
  ~SplitBackwardCorrectionData();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitBackwardCorrectionData(
      const SplitBackwardCorrectionData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitBackwardCorrectionData& operator=(
      const SplitBackwardCorrectionData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitBackwardCorrectionData(
      SplitBackwardCorrectionData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitBackwardCorrectionData& operator=(
      SplitBackwardCorrectionData&&) noexcept = default;

  void setImpulseStatus(const int dimi);

  Eigen::Block<Eigen::MatrixXd> KKT_mat_inv();

  const Eigen::Block<const Eigen::MatrixXd> KKT_mat_inv() const;

  Eigen::Block<Eigen::MatrixXd> auxMat();

  const Eigen::Block<const Eigen::MatrixXd> auxMat() const;

  Eigen::VectorBlock<Eigen::VectorXd> splitDirection();

  const Eigen::VectorBlock<const Eigen::VectorXd> splitDirection() const;

  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  Eigen::VectorBlock<Eigen::VectorXd> dxi();

  const Eigen::VectorBlock<const Eigen::VectorXd> dxi() const;

  Eigen::VectorBlock<Eigen::VectorXd> du();

  const Eigen::VectorBlock<const Eigen::VectorXd> du() const;

  Eigen::VectorBlock<Eigen::VectorXd> dq();

  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  Eigen::VectorBlock<Eigen::VectorXd> dv();

  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

private:
  int dimv_, dimx_, dimu_, dimi_, dimKKT_, du_begin_, dq_begin_, dv_begin_;
  Eigen::MatrixXd KKT_mat_inv_full_;
  Eigen::VectorXd split_direction_full_;

};

} // namespace idocp

#include "idocp/ocp/split_backward_correction_data.hxx"

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HPP_ 