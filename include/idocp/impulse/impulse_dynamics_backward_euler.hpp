#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler_data.hpp"


namespace idocp {

///
/// @class ImpulseDynamicsBackwardEuler
/// @brief Impulse dynamics constraint for the backward Euler.
///
class ImpulseDynamicsBackwardEuler {
public:
  ///
  /// @brief Constructs the impulse dynamics.
  /// @param[in] robot Robot model. 
  ///
  ImpulseDynamicsBackwardEuler(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseDynamicsBackwardEuler();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseDynamicsBackwardEuler();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseDynamicsBackwardEuler(
      const ImpulseDynamicsBackwardEuler&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseDynamicsBackwardEuler& operator=(
      const ImpulseDynamicsBackwardEuler&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseDynamicsBackwardEuler(
      ImpulseDynamicsBackwardEuler&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseDynamicsBackwardEuler& operator=(
      ImpulseDynamicsBackwardEuler&&) noexcept = default;

  ///
  /// @brief Linearizes the impulse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void linearizeImpulseDynamics(Robot& robot, 
                                const ImpulseStatus& impulse_status,  
                                const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the inverse impulse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] data Data for impulse dynamics.
  ///
  static void linearizeInverseImpulseDynamics(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const ImpulseSplitSolution& s, ImpulseDynamicsBackwardEulerData& data);

  ///
  /// @brief Linearizes the impulse velocity constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  static void linearizeImpulseVelocityConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual);

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
  /// @brief Condenses the inverse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in, out] data Data for impulse dynamics.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  static void condensing(const Robot& robot, 
                         ImpulseDynamicsBackwardEulerData& data, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this impulse stage.
  /// @param[in] robot Robot model. 
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void computeCondensedPrimalDirection(const Robot& robot, 
                                       const ImpulseSplitKKTMatrix& kkt_matrix, 
                                       ImpulseSplitDirection& d) const;

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this impulse stage.
  /// @param[in] robot Robot model. 
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void computeCondensedDualDirection(const Robot& robot, 
                                     ImpulseSplitDirection& d);

  ///
  /// @brief Computes the residual in the impulse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns l1-norm of the residual in the impulse dynamics constraint. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage.
  /// @return l1-norm of the residual in the impulse dynamics constraint.
  ///
  double l1NormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Returns squared norm of the residual in the impulse dynamics 
  /// constraint. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage.
  /// @return Squared norm of the residual in the impulse dynamics constraint.
  ///
  double squaredNormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ImpulseDynamicsBackwardEulerData data_;
  int dimv_, dimf_;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_backward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_ 