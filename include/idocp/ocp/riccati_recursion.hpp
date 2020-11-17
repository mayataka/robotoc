#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/hybrid/hybrid_container.hpp"

namespace idocp {

///
/// @class RiccatiRecursion 
/// @brief Contact sequence, i.e., sequence of contact status over the 
/// horizon. Provides the formulation of the optimal control problem 
/// with impulse and lift.
///
class RiccatiRecursion {
public:
  using HybridSplitOCPs = hybrid_container<SplitOCP, SplitImpulseOCP>;
  using HybridSplitSolutions = hybrid_container<SplitSolution, ImpulseSplitSolution>;
  using HybridSplitDirections = hybrid_container<SplitDirection, ImpulseSplitDirection>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  ///
  RiccatiRecursion(const Robot& robot, const double T, const int N, 
                   const int max_num_impulse=0, const int nproc=1);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiRecursion();

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiRecursion(const RiccatiRecursion&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  RiccatiRecursion& operator=(const RiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiRecursion(RiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiRecursion& operator=(RiccatiRecursion&&) noexcept = default;

  void backwardRiccatiRecursionTerminal(const TerminalOCP& terminal_ocp);

  void backwardRiccatiRecursion(const ContactSequence& contact_sequence, 
                                HybridSplitOCPs& split_ocps);

  void forwardRiccatiRecursionParallel(const ContactSequence& contact_sequence, 
                                       HybridSplitOCPs& split_ocps) const;

  void forwardRiccatiRecursionSerial(const ContactSequence& contact_sequence,
                                     HybridSplitOCPs& split_ocps);

  template <typename VectorType>
  void computeLagrangeMultiplierDirection(
      const ContactSequence& contact_sequence, const HybridSplitOCPs& split_ocps,
      const Eigen::MatrixBase<VectorType>& dx0, HybridSplitDirections& d);

  void computeDirectionAndStepSize(std::vector<Robot>& robots,
                                   const ContactSequence& contact_sequence,
                                   const HybridSplitOCPs& split_ocps,
                                   const HybridSplitSolutions& s,
                                   HybridSplitDirections& d);
  
  double primalStepSize() const;

  double dualStepSize() const;

private:
  int N_, max_num_impulse_, nproc_;
  double dtau_;
  std::vector<RiccatiFactorization> riccati_factorization_, 
                                    impulse_riccati_factorization_, 
                                    aux_riccati_factorization_,
                                    lift_riccati_factorization_;
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization_;
  StateConstraintRiccatiFactorizer constraint_factorizer_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_;

};

} // namespace idocp 

#include "idocp/ocp/riccati_recursion.hxx"

#endif // IDOCP_RICCATI_RECURSION_HPP_ 