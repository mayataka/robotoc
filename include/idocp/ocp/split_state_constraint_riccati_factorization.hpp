#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 
#define IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

///
/// @class SplitStateConstraintRiccatiFactorization
/// @brief Riccati factorized matrix and vector for a state only equality 
/// constraint split into each impulse stages.
///
class SplitStateConstraintRiccatiFactorization {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of possible impulse. Must be 
  /// more than 1 and less than N. 
  ///
  SplitStateConstraintRiccatiFactorization(const Robot& robot, const int N, 
                                           const int max_num_impulse);

  ///
  /// @brief Default constructor. 
  ///
  SplitStateConstraintRiccatiFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~SplitStateConstraintRiccatiFactorization();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitStateConstraintRiccatiFactorization(
      const SplitStateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitStateConstraintRiccatiFactorization& operator=(
      const SplitStateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitStateConstraintRiccatiFactorization(
      SplitStateConstraintRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitStateConstraintRiccatiFactorization& operator=(
      SplitStateConstraintRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Set the dimension of the impulse. 
  /// @param[in] impulse_status Impulse status. 
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set the dimension of the impulse 0. 
  ///
  void resetImpulseStatus();

  ///
  /// @brief A factorization matrix at a time stage. 
  /// @param[in] time_stage Time stage of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T(const int time_stage);

  ///
  /// @brief const version of SplitStateConstraintRiccatiFactorization::T(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T(const int time_stage) const;

  ///
  /// @brief A factorization matrix at a impulse time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_impulse(const int impulse_index);

  ///
  /// @brief const version of 
  /// SplitStateConstraintRiccatiFactorization::T_impulse(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_impulse(
      const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a aux time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_aux(const int impulse_index);

  ///
  /// @brief const version of 
  /// SplitStateConstraintRiccatiFactorization::T_aux(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_aux(const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a lift time stage. 
  /// @param[in] lift_index Lift index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_lift(const int lift_index);

  ///
  /// @brief const version of 
  /// SplitStateConstraintRiccatiFactorization::T_lift(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_lift(const int lift_index) const;

  ///
  /// @brief Partial derivative of the equality constriant with respect to the 
  /// configuration. 
  ///
  Eigen::Block<Eigen::MatrixXd> Eq();

  ///
  /// @brief const version of SplitStateConstraintRiccatiFactorization::Eq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> Eq() const;

  ///
  /// @brief Product of the partial derivative of the equality constriant with 
  /// respect to state and SplitRiccatiFactorization::N. 
  ///
  Eigen::Block<Eigen::MatrixXd> EN();

  ///
  /// @brief const version of SplitStateConstraintRiccatiFactorization::EN(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> EN() const;

  ///
  /// @brief Top Robot::dimv() rows of SplitRiccatiFactorization::EN(). 
  ///
  Eigen::Block<Eigen::MatrixXd> ENq();

  ///
  /// @brief const version of SplitStateConstraintRiccatiFactorization::ENq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENq() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::MatrixXd> T_full_, T_impulse_full_, T_aux_full_, T_lift_full_;
  Eigen::MatrixXd E_full_, EN_full_;
  int N_, max_num_impulse_, dimv_, dimx_, dimf_;

};

} // namespace idocp

#include "idocp/ocp/split_state_constraint_riccati_factorization.hxx"

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 