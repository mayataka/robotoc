#ifndef IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_
#define IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
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
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  RiccatiDirectionCalculator(const int N, const int max_num_impulse, 
                             const int nthreads);

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
                                           const KKTMatrix& kkt_matrix, 
                                           const Solution& s, Direction& d);

  void computeNewtonDirectionFromRiccatiFactorization(
      OCP& ocp, std::vector<Robot>& robots, 
      const RiccatiFactorization& factorization, 
      const Solution& s, Direction& d);

  double maxPrimalStepSize() const;

  double maxDualStepSize() const;

private:
  int N_, nthreads_, N_all_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace idocp 

#endif // IDOCP_RICCATI_DIRECTION_CALCULATOR_HPP_ 