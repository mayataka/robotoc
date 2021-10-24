#ifndef ROBOTOC_IMPULSE_DYNAMICS_HPP_
#define ROBOTOC_IMPULSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_dynamics_data.hpp"


namespace robotoc {

///
/// @class ImpulseDynamics
/// @brief Impulse dynamics constraint for the forward Euler.
///
class ImpulseDynamics {
public:
  ///
  /// @brief Constructs the impulse dynamics.
  /// @param[in] robot Robot model. 
  ///
  ImpulseDynamics(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseDynamics();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseDynamics();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseDynamics(const ImpulseDynamics&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseDynamics& operator=(const ImpulseDynamics&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseDynamics(ImpulseDynamics&&) noexcept = default;
  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseDynamics& operator=(ImpulseDynamics&&) noexcept = default;

  ///
  /// @brief Computes the residual in the impulse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  ///
  void evalImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status,
                           const ImpulseSplitSolution& s);

  ///
  /// @brief Linearizes the impulse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void linearizeImpulseDynamics(Robot& robot, 
                                const ImpulseStatus& impulse_status, 
                                const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Condenses the inverse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void condenseImpulseDynamics(Robot& robot, 
                               const ImpulseStatus& impulse_status,
                               ImpulseSplitKKTMatrix& kkt_matrix, 
                               ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Expands the primal variables, i.e., computes the Newton direction 
  /// of the condensed primal variables (impulse change in the velocity dv and 
  /// the impulse forces f) of this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void expandPrimal(ImpulseSplitDirection& d) const;

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variables (Lagrange multipliers) of this impulse 
  /// stage.
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  template <typename SplitDirectionType>
  void expandDual(const SplitDirectionType& d_next, ImpulseSplitDirection& d);

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual of the impulse dynamics constraint. 
  /// @return Squared norm of the KKT residual in the impulse dynamics 
  /// constraint.
  ///
  double KKTError() const;

  ///
  /// @brief Returns l1-norm of the constraint violation, that is, the primal
  /// residual in the impulse dynamics constraint. 
  /// @return l1-norm of the constraint violation.
  ///
  double constraintViolation() const;

private:
  ImpulseDynamicsData data_;

};

} // namespace robotoc 

#include "robotoc/impulse/impulse_dynamics.hxx"

#endif // ROBOTOC_IMPULSE_DYNAMICS_HPP_ 