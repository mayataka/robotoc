#ifndef IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HPP_
#define IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/contact_sequence.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"

namespace idocp {

class StateConstraintsRiccatiFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  StateConstraintsRiccatiFactorizer(const Robot& robot, const int N, 
                                    const int nproc);

  // Default constructor.
  StateConstraintsRiccatiFactorizer();

  // Destructor.
  ~StateConstraintsRiccatiFactorizer();
 
  // Default copy constructor.
  StateConstraintsRiccatiFactorizer(
      const StateConstraintsRiccatiFactorizer&) = default;

  // Default copy operator.
  StateConstraintsRiccatiFactorizer& operator=(
      const StateConstraintsRiccatiFactorizer&) = default;

  // Default move constructor.
  StateConstraintsRiccatiFactorizer(
      StateConstraintsRiccatiFactorizer&&) noexcept = default;

  // Default move assign operator.
  StateConstraintsRiccatiFactorizer& operator=(
      StateConstraintsRiccatiFactorizer&&) noexcept = default;

  void computeLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence,
      const std::vector<StateConstraintRiccatiFactorization>& factorizations,
      const std::vector<RiccatiSolution>& riccati_solutions,
      std::vector<ImpulseSplitDirection>& d);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::LLT<Eigen::MatrixXd>> llts_;
  int N_, dimv_, dimx_, nproc_;

};

} // namespace idocp

#include "idocp/ocp/state_constraints_riccati_factorizer.hxx"

#endif // IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZATER_HPP_ 