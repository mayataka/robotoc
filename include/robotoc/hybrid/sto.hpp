#ifndef ROBOTOC_STO_HPP_
#define ROBOTOC_STO_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"


namespace robotoc {

///
/// @class STO 
/// @brief The switching time optimization (STO) problem.
///
class STO {
public:
  ///
  /// @brief Construct the STO problem. 
  /// @param[in] sto_cost Shared ptr to the STO cost function.
  /// @param[in] sto_constraints Shared ptr to the STO constraints.
  /// @param[in] max_num_impulse_events Maximum number of the impulse on the 
  /// horizon. Must be non-negative. 
  ///
  STO(const std::shared_ptr<STOCostFunction>& sto_cost, 
      const std::shared_ptr<STOConstraints>& sto_constraints,
      const int max_num_impulse_events);

  ///
  /// @brief Default Constructor.
  ///
  STO();

  ///
  /// @brief Default copy constructor. 
  ///
  STO(const STO&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  STO& operator=(const STO&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STO(STO&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STO& operator=(STO&&) noexcept = default;

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] ocp Optimal control problem.
  ///
  void initConstraints(const OCP& ocp) const;

  ///
  /// @brief Computes the KKT residual of switching time otpimization problem. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTResidual(const OCP& ocp, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT system, i.e., the condensed KKT matrix and KKT
  /// residual for Newton's method. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTSystem(const OCP& ocp, KKTMatrix& kkt_matrix, 
                        KKTResidual& kkt_residual);

  double KKTError() const;

  double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

  double totalCost() const;

  ///
  /// @brief Integrates the solution (switching times). The internal 
  /// switching times of contact_sequence are updated. The slack and dual
  /// variables in STO constraints are also updated.
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  ///
  void integrateSolution(const OCP& ocp, 
                         std::shared_ptr<ContactSequence>& contact_sequence,
                         const double primal_step_size,
                         const double dual_step_size,
                         const Direction& d) const;

private:
  std::shared_ptr<STOCostFunction> sto_cost_;
  std::shared_ptr<STOConstraints> sto_constraints_;
  int max_num_impulse_events_;
  double kkt_error_, cost_val_;
  Eigen::VectorXd h_phase_;
  bool is_sto_enabled_;
};

} // namespace robotoc

#include "robotoc/hybrid/sto.hxx"

#endif // ROBOTOC_STO_HPP_ 