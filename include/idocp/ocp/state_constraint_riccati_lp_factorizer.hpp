#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HPP_ 
#define IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HPP_ 

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"

namespace idocp {

///
/// @class StateConstraintRiccatiLPFactorizer 
/// @brief A factorizer for Linear programming for the Riccati recursion with
/// pure-state constraints.
///
class StateConstraintRiccatiLPFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  StateConstraintRiccatiLPFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  StateConstraintRiccatiLPFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~StateConstraintRiccatiLPFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  StateConstraintRiccatiLPFactorizer(
      const StateConstraintRiccatiLPFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  StateConstraintRiccatiLPFactorizer& operator=(
      const StateConstraintRiccatiLPFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateConstraintRiccatiLPFactorizer(
      StateConstraintRiccatiLPFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateConstraintRiccatiLPFactorizer& operator=(
      StateConstraintRiccatiLPFactorizer&&) noexcept = default;

  ///
  /// @brief Factorize matrices and vectors for the linear problem to obtain
  /// the directions of the Lagrange multipliers. Used in
  /// StateConstraintRiccatiFactorizer::computeLagrangeMultiplierDirection().
  /// @param[in] impulse_riccati_factorization Riccati factorizations for 
  /// an impulse stage.
  /// @param[in, out] constraint_factorization A constraint factorization.
  /// @param[in] dx0 Newton direction at the initial stage.
  /// @param[in] constraint_index Constraint index.
  ///
  template <typename VectorType>
  void factorizeLinearProblem(
      const ContactSequence& constact_sequence,
      const RiccatiFactorization& impulse_riccati_factorization,
      StateConstraintRiccatiFactorization& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0,
      const int constraint_index);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_;
  Eigen::VectorXd Pidx0q_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HPP_ 