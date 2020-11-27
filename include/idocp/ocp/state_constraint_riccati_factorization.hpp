#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_state_constraint_riccati_factorization.hpp"
#include "idocp/hybrid/contact_sequence.hpp"

namespace idocp {

///
/// @class StateConstraintRiccatiFactorization
/// @brief Riccati factorized matrix and vector for the pure-state equality 
/// constraints.
///
class StateConstraintRiccatiFactorization {
public:
  ///
  /// @brief Allocate factorization matrices and vectors.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of possible impulse. Must be 
  /// more than 1 and less than N. 
  ///
  StateConstraintRiccatiFactorization(const Robot& robot, const int N, 
                                      const int max_num_impulse);

  ///
  /// @brief Default constructor. 
  ///
  StateConstraintRiccatiFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~StateConstraintRiccatiFactorization();
 
  ///
  /// @brief Default copy constructor. 
  ///
  StateConstraintRiccatiFactorization(
      const StateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  StateConstraintRiccatiFactorization& operator=(
      const StateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateConstraintRiccatiFactorization(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateConstraintRiccatiFactorization& operator=(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Set the dimension of the constraints. 
  /// @param[in] contact_sequence Contact sequence. 
  ///
  void setContactSequence(const ContactSequence& contact_sequence);

  ///
  /// @brief A factorization matrix at a time stage. 
  /// @param[in] time_stage Time stage of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T(const int constraint_index, 
                                  const int time_stage);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T(const int constraint_index, 
                                              const int time_stage) const;

  ///
  /// @brief A factorization matrix at a impulse time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_impulse(const int constraint_index,
                                          const int impulse_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_impulse(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_impulse(
      const int constraint_index, const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a aux time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_aux(const int constraint_index, 
                                      const int impulse_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_aux(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_aux(
      const int constraint_index, const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a lift time stage. 
  /// @param[in] lift_index Lift index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_lift(const int constraint_index, 
                                       const int lift_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_lift(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_lift(const int constraint_index, 
                                                   const int lift_index) const;

  ///
  /// @brief Partial derivative of the equality constriant with respect to the 
  /// configuration. 
  ///
  Eigen::Block<Eigen::MatrixXd> Eq(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::Eq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> Eq(
      const int constraint_index) const;

  ///
  /// @brief Product of the partial derivative of the equality constriant with 
  /// respect to state and RiccatiFactorization::N. 
  ///
  Eigen::Block<Eigen::MatrixXd> EN(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::EN(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> EN(
      const int constraint_index) const;

  ///
  /// @brief Top Robot::dimv() rows of RiccatiFactorization::EN(). 
  ///
  Eigen::Block<Eigen::MatrixXd> ENq(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENq(
      const int constraint_index) const;

  ///
  /// @brief Product of RiccatiFactorization::EN() and 
  /// StateConstraintRiccatiFactorization::Eq(). 
  /// @param[in] constraint_index Constraint index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> ENEt(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENEt(). 
  /// @param[in] constraint_index Constraint index of interested. 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENEt(
      const int constraint_index) const;

  ///
  /// @brief Product of RiccatiFactorization::EN() and 
  /// StateConstraintRiccatiFactorization::Eq(). 
  /// @param[in] constraint_index Constraint index of interested. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> ENT(const int constraint_index, 
                                    const int impulse_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENEt(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENT(
      const int constraint_index, const int impulse_index) const;

  ///
  /// @brief Product of RiccatiFactorization::EN() and 
  /// StateConstraintRiccatiFactorization::Eq(). 
  /// @param[in] constraint_index Constraint index of interested. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> ENT();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENT(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENT() const;

  ///
  /// @brief Residual of the equality constriant.
  /// @param[in] constraint_index Constraint index of interested. 
  ///
  Eigen::VectorBlock<Eigen::VectorXd> e(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::e(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> e(
      const int constraint_index) const;

  ///
  /// @brief Residual of the equality constriant.
  /// @param[in] constraint_index Constraint index of interested. 
  ///
  Eigen::VectorBlock<Eigen::VectorXd> e();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::e(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> e() const;

  ///
  /// @brief Direction of the pure-state equality constraint.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dxi();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::dxi(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dxi() const;

  ///
  /// @brief Direction of the pure-state equality constraint.
  /// @param[in] constraint_index Constraint index of interested. 
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dxi(const int constraint_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::dxi(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dxi(
      const int constraint_index) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd ENT_full_;
  Eigen::VectorXd e_full_, dxi_full_;
  std::vector<SplitStateConstraintRiccatiFactorization> factorizations_;
  std::vector<int> dimf_, f_begin_;
  int N_, max_num_impulse_, dimv_, dimx_, dimf_total_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorization.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 