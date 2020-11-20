#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPLinearizer
/// @brief Linearize of the optimal control problem. 
///
class OCPLinearizer {
public:
  using HybridOCP = hybrid_container<SplitOCP, SplitImpulseOCP>;
  using HybridSolution = hybrid_container<SplitSolution, ImpulseSplitSolution>;
  using HybridKKTMatrix = hybrid_container<KKTMatrix, ImpulseKKTMatrix>;
  using HybridKKTResidual = hybrid_container<KKTResidual, ImpulseKKTResidual>;

  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPLinearizer(const double T, const int N, const int max_num_impulse=0, 
                const int num_proc=1);

  ///
  /// @brief Default constructor. 
  ///
  OCPLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPLinearizer(const OCPLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPLinearizer& operator=(const OCPLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPLinearizer(OCPLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPLinearizer& operator=(OCPLinearizer&&) noexcept = default;

  void linearizeOCP(HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
                    std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, const HybridSolution& s,
                    HybridKKTMatrix& kkt_matrix,
                    HybridKKTResidual& kkt_residual);

  void computeKKTResidual(HybridOCP& split_ocps, TerminalOCP& terminal_ocp,
                          std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const HybridSolution& s,
                          HybridKKTMatrix& kkt_matrix, 
                          HybridKKTResidual& kkt_residual);

  // void computeDirections(
  //     std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
  //     const HybridRiccatiFactorizer& factorizer, 
  //     HybridRiccatiFactorization& factorization, 
  //     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization, 
  //     const HybridSolution& s, HybridDirection& d);

  // double maxPrimalStepSize(const ContactSequence& contact_sequence) const;

  // double maxDualStepSize(const ContactSequence& contact_sequence) const;

  // void updateSolution(const std::vector<Robot>& robots, 
  //                     const HybridKKTMatrix& kkt_matrix,
  //                     const HybridKKTResidual& kkt_residual,
  //                     const double primal_step_size, 
  //                     const double dual_step_size, HybridDirection& d, 
  //                     HybridSolution& s);

  // static void computePrimalDirectionInitial(
  //     const RiccatiFactorizer factorizer, 
  //     const RiccatiFactorization factorization, SplitDirection& d, 
  //     const bool exist_state_constraint);

  // template <typename VectorType>
  // static void computePrimalDirection(const RiccatiFactorizer factorizer,
  //                                    const RiccatiFactorization factorization,
  //                                    const Eigen::MatrixBase<VectorType>& dx0, 
  //                                    SplitDirection& d,
  //                                    const bool exist_state_constraint);

  // template <typename VectorType>
  // static void computePrimalDirectionTerminal(
  //     const RiccatiFactorization factorization, 
  //     const Eigen::MatrixBase<VectorType>& dx0, SplitDirection& d);

  // template <typename VectorType>
  // static void computePrimalDirectionImpulse(
  //     const RiccatiFactorization factorization, 
  //     const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d);

  // static void aggregateLagrangeMultiplierDirection(
  //     const ContactSequence& contact_sequence,
  //     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
  //     const std::vector<ImpulseSplitDirection>& d_impulse, const int time_stage,
  //     RiccatiFactorization& riccati_factorization);

  // static void aggregateLagrangeMultiplierDirectionImpulse(
  //     const ContactSequence& contact_sequence,
  //     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
  //     const std::vector<ImpulseSplitDirection>& d_impulse, 
  //     const int impulse_index,
  //     RiccatiFactorization& impulse_riccati_factorization);

  // static void aggregateLagrangeMultiplierDirectionAux(
  //     const ContactSequence& contact_sequence,
  //     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
  //     const std::vector<ImpulseSplitDirection>& d_impulse, 
  //     const int impulse_index,
  //     RiccatiFactorization& aux_riccati_factorization);

  // static void aggregateLagrangeMultiplierDirectionLift(
  //     const ContactSequence& contact_sequence,
  //     const std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
  //     const std::vector<ImpulseSplitDirection>& d_impulse, 
  //     const int lift_index,
  //     RiccatiFactorization& lift_riccati_factorization);

private:

  template <typename Algorithm>
  void runParallel(HybridOCP& split_ocps, TerminalOCP& terminal_ocp,
                   std::vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const HybridSolution& s,
                   HybridKKTMatrix& kkt_matrix, HybridKKTResidual& kkt_residual);

  const Eigen::VectorXd& q_prev(const ContactSequence& contact_sequence, 
                                const Eigen::VectorXd& q,
                                const HybridSolution& s,
                                const int time_stage) const;

  double dtau(const ContactSequence& contact_sequence, 
              const int time_stage) const;

  double T_, dtau_;
  int N_, num_proc_;

};

} // namespace idocp 

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_