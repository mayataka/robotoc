#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"
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
                                   const int nproc);

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
  /// @param[in] contact_sequence Contact sequence.
  /// @param[in] impulse_riccati_factorization Riccati factorizations for 
  /// impulse stages.
  /// @param[in, out] constraint_factorization Constraint factorizations.
  /// @param[in] dx0 Newton direction at the initial stage.
  /// @param[in, out] d_impulse Directions at the impulse stages.
  ///
  template <typename VectorType>
  void computeLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence,
      const std::vector<RiccatiFactorization>& impulse_riccati_factorization,
      std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0,
      std::vector<ImpulseSplitDirection>& d_impulse);


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::LLT<Eigen::MatrixXd>> llts_;
  std::vector<StateConstraintRiccatiLPFactorizer> lp_factorizer_;
  int max_num_impulse_, dimv_, nproc_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorizer.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATER_HPP_ 