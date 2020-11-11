#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/contact_sequence.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"

namespace idocp {

class StateConstraintRiccatiFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  StateConstraintRiccatiFactorizer(const Robot& robot, 
                                   const int max_num_impulse, const int nproc);

  // Default constructor.
  StateConstraintRiccatiFactorizer();

  // Destructor.
  ~StateConstraintRiccatiFactorizer();
 
  // Default copy constructor.
  StateConstraintRiccatiFactorizer(
      const StateConstraintRiccatiFactorizer&) = default;

  // Default copy operator.
  StateConstraintRiccatiFactorizer& operator=(
      const StateConstraintRiccatiFactorizer&) = default;

  // Default move constructor.
  StateConstraintRiccatiFactorizer(
      StateConstraintRiccatiFactorizer&&) noexcept = default;

  // Default move assign operator.
  StateConstraintRiccatiFactorizer& operator=(
      StateConstraintRiccatiFactorizer&&) noexcept = default;

  template <typename VectorType>
  void computeLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence,
      const std::vector<RiccatiFactorization>& impulse_riccati_factorization,
      std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0,
      std::vector<ImpulseSplitDirection>& d);

  template <typename VectorType>
  void factorizeLinearProblem(
      const RiccatiFactorization& impulse_riccati_factorization,
      StateConstraintRiccatiFactorization& constraint_factorization,
      const Eigen::MatrixBase<VectorType>& dx0);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::LLT<Eigen::MatrixXd>> llts_;
  int max_num_impulse_, dimv_, nproc_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorizer.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATER_HPP_ 