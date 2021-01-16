#ifndef IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 
#define IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class PreImpulseSplitKKTResidual
/// @brief KKT residual split at the impulse stage. 
///
class PreImpulseSplitKKTResidual {
public:
  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  PreImpulseSplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  PreImpulseSplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~PreImpulseSplitKKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  PreImpulseSplitKKTResidual(const PreImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  PreImpulseSplitKKTResidual& operator=(
      const PreImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  PreImpulseSplitKKTResidual(PreImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  PreImpulseSplitKKTResidual& operator=(
      PreImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Set impulse status from robot model, i.e., set dimension of the 
  /// impulse and equality constraints.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief Set the KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Chech the equivalence of two PreImpulseSplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const PreImpulseSplitKKTResidual& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd P_full_;
  int dimf_;

};

} // namespace idocp 

#include "idocp/ocp/pre_impulse_split_kkt_residual.hxx"

#endif // IDOCP_PRE_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 