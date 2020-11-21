#ifndef IDOCP_OCP_SOLUTION_INTEGRATOR_HPP_
#define IDOCP_OCP_SOLUTION_INTEGRATOR_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {

///
/// @class OCPSolutionIntegrator
/// @brief Linearize of the optimal control problem. 
///
class OCPSolutionIntegrator {
public:
  using HybridOCP = hybrid_container<SplitOCP, SplitImpulseOCP>;
  using HybridSolution = hybrid_container<SplitSolution, ImpulseSplitSolution>;
  using HybridDirection = hybrid_container<SplitDirection, ImpulseSplitDirection>;
  using HybridKKTMatrix = hybrid_container<KKTMatrix, ImpulseKKTMatrix>;
  using HybridKKTResidual = hybrid_container<KKTResidual, ImpulseKKTResidual>;

  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPSolutionIntegrator(const double T, const int N, 
                        const int max_num_impulse=0, const int num_proc=1);

  ///
  /// @brief Default constructor. 
  ///
  OCPSolutionIntegrator();

  ///
  /// @brief Destructor. 
  ///
  ~OCPSolutionIntegrator();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPSolutionIntegrator(const OCPSolutionIntegrator&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPSolutionIntegrator& operator=(const OCPSolutionIntegrator&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPSolutionIntegrator(OCPSolutionIntegrator&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPSolutionIntegrator& operator=(OCPSolutionIntegrator&&) noexcept = default;

  void integrate(HybridOCP& split_ocps, TerminalOCP& terminal_ocp, 
                 std::vector<Robot>& robots,
                 const ContactSequence& contact_sequence,
                 const HybridKKTMatrix& kkt_matrix,
                 const HybridKKTResidual& kkt_residual,
                 const double primal_step_size,
                 const double dual_step_size,
                 HybridDirection& d, HybridSolution& s) const;

private:

  double T_, dtau_;
  int N_, num_proc_;

};

} // namespace idocp 

#include "idocp/ocp/ocp_solution_integrator.hxx"

#endif // IDOCP_OCP_SOLUTION_INTEGRATOR_HPP_ 