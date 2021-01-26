#ifndef IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HPP_
#define IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class ImpulseSplitBackwardCorrectionData
/// @brief Split unconstrained backward correction.
///
class ImpulseSplitBackwardCorrectionData {
public:
  ///
  /// @brief Constructor.
  ///
  ImpulseSplitBackwardCorrectionData(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitBackwardCorrectionData();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitBackwardCorrectionData();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitBackwardCorrectionData(
      const ImpulseSplitBackwardCorrectionData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitBackwardCorrectionData& operator=(
      const ImpulseSplitBackwardCorrectionData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitBackwardCorrectionData(
      ImpulseSplitBackwardCorrectionData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitBackwardCorrectionData& operator=(
      ImpulseSplitBackwardCorrectionData&&) noexcept = default;

  void setImpulseStatus(const int dimf);

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

  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

  Eigen::VectorBlock<Eigen::VectorXd> df();

  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  Eigen::VectorBlock<Eigen::VectorXd> dq();

  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  Eigen::VectorBlock<Eigen::VectorXd> dv();

  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

private:
  int dimv_, dimx_, dimf_, dimKKT_, df_begin_, dq_begin_, dv_begin_;
  Eigen::MatrixXd KKT_mat_inv_full_;
  Eigen::VectorXd split_direction_full_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_backward_correction_data.hxx"

#endif // IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HPP_ 