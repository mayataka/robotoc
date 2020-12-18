#ifndef IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_
#define IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/split_riccati_factorizer.hpp"
#include "idocp/impulse/impulse_split_riccati_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"


namespace idocp {

///
/// @class RiccatiDirectionCalculator 
/// @brief Linearize of the optimal control problem. 
///
class RiccatiDirectionCalculator {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  RiccatiDirectionCalculator(const int N, const int max_num_impulse, 
                             const int num_proc);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiDirectionCalculator();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiDirectionCalculator();

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiDirectionCalculator(const RiccatiDirectionCalculator&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  RiccatiDirectionCalculator& operator=(
      const RiccatiDirectionCalculator&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiDirectionCalculator(RiccatiDirectionCalculator&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiDirectionCalculator& operator=(
      RiccatiDirectionCalculator&&) noexcept = default;

  static void computeInitialStateDirection(const std::vector<Robot>& robots, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v, 
                                           const Solution& s, Direction& d);

  void computeNewtonDirectionFromRiccatiFactorization(
      OCP& ocp, const OCPDiscretizer& ocp_discretizer,
      std::vector<Robot>& robots, const RiccatiFactorizer& factorizer, 
      const RiccatiFactorization& factorization, 
      const Solution& s, Direction& d);

  double maxPrimalStepSize() const;

  double maxDualStepSize() const;

  static const SplitRiccatiFactorization& next_riccati_factorization(
      const OCPDiscretizer& ocp_discretizer, 
      const RiccatiFactorization& factorization, const int time_stage) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(time_stage)) {
      return factorization.impulse[ocp_discretizer.impulseIndex(time_stage)];
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(time_stage)) {
      return factorization.lift[ocp_discretizer.liftIndex(time_stage)];
    }
    else {
      return factorization[time_stage+1];
    }
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int N_, num_proc_, N_all_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace idocp 

#endif // IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_ 