#ifndef IDOCP_UNRICCATI_SOLVER_HPP_
#define IDOCP_UNRICCATI_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/unconstrained_container.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"
#include "idocp/unocp/unconstrained_container.hpp"

namespace idocp {

///
/// @class RiccatiSolver
/// @brief Riccati recursion.
///
class UnRiccatiSolver {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nproc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  UnRiccatiSolver(const Robot& robot, const double T, const int N, 
                  const int nproc);

  ///
  /// @brief Default constructor. 
  ///
  UnRiccatiSolver();

  ///
  /// @brief Destructor. 
  ///
  ~UnRiccatiSolver();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnRiccatiSolver(const UnRiccatiSolver&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnRiccatiSolver& operator=(const UnRiccatiSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnRiccatiSolver(UnRiccatiSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnRiccatiSolver& operator=(UnRiccatiSolver&&) noexcept = default;

  void linearizeOCP(UnOCP& ocp, const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v);

  ///
  /// @brief Performs the backward Riccati recursion for the terminal stage. 
  /// @param[in] terminal_kkt_matrix KKT matrix at the terminal stage. 
  /// @param[in] terminal_kkt_residual KKT residual at the terminal stage. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiSolverTerminal(
      const SplitKKTMatrix& terminal_kkt_matrix, 
      const SplitKKTResidual& terminal_kkt_residual,
      UnRiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the backward Riccati recursion. Call 
  /// RiccatiSolver::backwardRiccatiSolverTerminal() before calling this
  /// function.
  /// @param[in, out] unkkt_matrix KKT matrix. 
  /// @param[in, out] unkkt_residual KKT residual. 
  ///
  void backwardRiccatiSolver(UnKKTMatrix& unkkt_matrix, 
                                UnKKTResidual& unkkt_residual,
                                UnRiccatiFactorization& ricccati_factorization);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] unkkt_residual KKT residual. 
  /// @param[in, out] d Direction. d[0].dx() must be computed before 
  /// calling this function.
  ///
  void forwardRiccatiSolver(const UnKKTResidual& unkkt_residual, 
                               UnDirection& d) const;

private:
  int N_;
  double T_, dtau_;
  UnRiccatiFactorizer factorizer_;

};

} // namespace idocp

#endif // IDOCP_UNRICCATI_SOLVER_HPP_ 