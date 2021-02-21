#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"

namespace idocp {

///
/// @class RiccatiRecursion
/// @brief Riccati recursion.
///
class RiccatiRecursion {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  ///
  RiccatiRecursion(const Robot& robot, const int N, const int max_num_impulse);

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
  /// @brief Default copy operator. 
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

  ///
  /// @brief Performs the backward Riccati recursion. Call 
  /// RiccatiRecursion::backwardRiccatiRecursionTerminal() before calling this
  /// function.
  /// @param[in, out] factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(const OCPDiscretizer& ocp_discretizer, 
                                KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
                                const StateConstraintJacobian& jac,
                                RiccatiFactorization& factorization);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Split direction. d[0].dx() must be computed before 
  /// calling this function.
  ///
  void forwardRiccatiRecursion(const OCPDiscretizer& ocp_discretizer, 
                               const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual, 
                               Direction& d) const;

private:
  RiccatiFactorizer factorizer_;

};

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HPP_ 