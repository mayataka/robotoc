#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"

namespace idocp {

///
/// @class StateConstraintRiccatiFactorizer 
/// @brief State constriant factorizer for Riccati recursion.
///
class StateConstraintRiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] max_num_impulse_phase Maximum number of impulse phases over the horizon.
  /// @param[in] nproc Number of threads used in this class.
  ///
  StateConstraintRiccatiFactorizer(const Robot& robot, 
                                   const int max_num_impulse_phase, 
                                   const int nproc=1);

  ///
  /// @brief Default constructor. 
  ///
  StateConstraintRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~StateConstraintRiccatiFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  StateConstraintRiccatiFactorizer(
      const StateConstraintRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  StateConstraintRiccatiFactorizer& operator=(
      const StateConstraintRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateConstraintRiccatiFactorizer(
      StateConstraintRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateConstraintRiccatiFactorizer& operator=(
      StateConstraintRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Computes the directions of the Lagrange multipliers with respect
  /// to the pure-state constraints.
  ///
  template <typename VectorType>
  void computeLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence,
      const std::vector<RiccatiFactorization>& impulse_riccati_factorization,
      std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0,
      std::vector<ImpulseSplitDirection>& d_impulse);

  ///
  /// @brief Factorize matrices and vectors for the linear problem to obtain
  /// the directions of the Lagrange multipliers.
  ///
  template <typename VectorType>
  void factorizeLinearProblem(
      const RiccatiFactorization& impulse_riccati_factorization,
      StateConstraintRiccatiFactorization& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0) const;

  void aggregateLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence,
      const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const std::vector<ImpulseSplitDirection>& d_impulse, const int time_stage,
      RiccatiFactorization& riccati_factorization) const;

  void aggregateLagrangeMultiplierDirectionImpulse(
      const ContactSequence& contact_sequence,
      const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const std::vector<ImpulseSplitDirection>& d_impulse, 
      const int impulse_index,
      RiccatiFactorization& riccati_factorization) const;

  void aggregateLagrangeMultiplierDirectionLift(
      const ContactSequence& contact_sequence,
      const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const std::vector<ImpulseSplitDirection>& d_impulse, 
      const int lift_index,
      RiccatiFactorization& riccati_factorization) const;


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::LLT<Eigen::MatrixXd>> llts_;
  int max_num_impulse_, dimv_, nproc_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorizer.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATER_HPP_ 