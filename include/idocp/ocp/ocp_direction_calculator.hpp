#ifndef IDOCP_DIRECTION_CALCULATOR_HPP_
#define IDOCP_DIRECTION_CALCULATOR_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPDirectionCalculator
/// @brief Linearize of the optimal control problem. 
///
class OCPDirectionCalculator {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPDirectionCalculator(const double T, const int N, 
                         const int max_num_impulse, const int num_proc);

  ///
  /// @brief Default constructor. 
  ///
  OCPDirectionCalculator();

  ///
  /// @brief Destructor. 
  ///
  ~OCPDirectionCalculator();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPDirectionCalculator(const OCPDirectionCalculator&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPDirectionCalculator& operator=(const OCPDirectionCalculator&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPDirectionCalculator(OCPDirectionCalculator&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPDirectionCalculator& operator=(OCPDirectionCalculator&&) noexcept = default;

  static void computeInitialStateDirection(const std::vector<Robot>& robots, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v, 
                                           const HybridSolution& s, 
                                           HybridDirection& d);

  void computeDirection(HybridOCP& split_ocps, std::vector<Robot>& robots, 
                        const ContactSequence& contact_sequence, 
                        const HybridRiccatiFactorizer& factorizer, 
                        HybridRiccatiFactorization& factorization, 
                        const HybridSolution& s, HybridDirection& d);

  double maxPrimalStepSize(const ContactSequence& contact_sequence) const;

  double maxDualStepSize(const ContactSequence& contact_sequence) const;

  static void computePrimalDirection(
      const RiccatiFactorizer factorizer, 
      const RiccatiFactorization factorization, 
      const RiccatiFactorization factorization_next, SplitDirection& d, 
      const bool exist_state_constraint);

  static void computePrimalDirectionTerminal(
      const RiccatiFactorization factorization, SplitDirection& d);

  static void computePrimalDirectionImpulse(
      const RiccatiFactorization factorization, ImpulseSplitDirection& d,
      const bool exist_state_constraint);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  double dtau(const ContactSequence& contact_sequence, 
              const int time_stage) const;

  double T_, dtau_;
  int N_, num_proc_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace idocp 

#include "idocp/ocp/ocp_direction_calculator.hxx"

#endif // IDOCP_DIRECTION_CALCULATOR_HPP_ 