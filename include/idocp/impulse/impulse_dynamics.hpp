#ifndef IDOCP_IMPULSE_DYNAMICS_HPP_
#define IDOCP_IMPULSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_data.hpp"


namespace idocp {

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
      const ImpulseSplitSolution& s, ImpulseDynamicsData& data);

  ///
  /// @brief Linearizes the impulse velocity constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in, out] data Data for impulse dynamics.
  ///
  static void linearizeImpulseVelocityConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseDynamicsData& data);

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
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in, out] data Data for impulse dynamics.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  static void condensing(const Robot& robot, 
                         const ImpulseStatus& impulse_status,
                         ImpulseDynamicsData& data, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void computeCondensedPrimalDirection(ImpulseSplitDirection& d) const;

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this impulse stage.
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in] kkt_residual Split KKT residual of this impulse stage.
  /// @param[in] dgmm Direction of the costate of the next time stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  template <typename VectorType>
  void computeCondensedDualDirection(const ImpulseSplitKKTMatrix& kkt_matrix, 
                                     const ImpulseSplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this impulse stage.
  /// @param[in] data Data for impulse dynamics.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  static void expandPrimal(const ImpulseDynamicsData& data, 
                           ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this impulse stage.
  /// @param[in] data Data for impulse dynamics.
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in] kkt_residual Split KKT residual of this impulse stage.
  /// @param[in] dgmm Direction of the costate of the next time stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  template <typename VectorType>
  static void expansionDual(ImpulseDynamicsData& data,
                            const ImpulseSplitKKTMatrix& kkt_matrix, 
                            const ImpulseSplitKKTResidual& kkt_residual,
                            const Eigen::MatrixBase<VectorType>& dgmm,
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

private:
  ImpulseDynamicsData data_;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_HPP_ 