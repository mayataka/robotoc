#ifndef IDOCP_UNCONSTRAINED_DYNAMICS_HPP_
#define IDOCP_UNCONSTRAINED_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"


namespace idocp {

///
/// @class UnconstrainedDynamics
/// @brief Inverse dynamics constraint without constraints (floating base or 
/// contacts).
///
class UnconstrainedDynamics {
public:
  ///
  /// @brief Constructs the inverse dynamics.
  /// @param[in] robot Robot model. 
  ///
  UnconstrainedDynamics(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrainedDynamics();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrainedDynamics();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrainedDynamics(const UnconstrainedDynamics&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrainedDynamics& operator=(const UnconstrainedDynamics&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrainedDynamics(UnconstrainedDynamics&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrainedDynamics& operator=(UnconstrainedDynamics&&) noexcept = default;

  ///
  /// @brief Linearizes the unconstrained dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeUnconstrainedDynamics(Robot& robot, const double dt, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual);

  ///
  /// @brief Condenses the unconstrained dynamics constraint. 
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] unkkt_matrix Condensed split KKT matrix of this time 
  /// stage.
  /// @param[in, out] unkkt_residual Condensed split KKT residual of this time 
  /// stage.
  ///
  void condenseUnconstrainedDynamics(const SplitKKTMatrix& kkt_matrix, 
                                     const SplitKKTResidual& kkt_residual,
                                     SplitUnKKTMatrix& unkkt_matrix, 
                                     SplitUnKKTResidual& unkkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed variables of this 
  /// time stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedDirection(const double dt, 
                                 const SplitKKTMatrix& kkt_matrix, 
                                 const SplitKKTResidual& kkt_residual, 
                                 SplitDirection& d);

  ///
  /// @brief Computes the residual in the unconstrained dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this time stage.
  ///
  void computeUnconstrainedDynamicsResidual(Robot& robot, 
                                            const SplitSolution& s);

  ///
  /// @brief Returns l1-norm of the residual in the unconstrained dynamics 
  /// constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return l1-norm of the residual in the unconstrained dynamics constraint.
  ///
  double l1NormUnconstrainedDynamicsResidual(const double dt) const;

  ///
  /// @brief Returns squared norm of the residual in the unconstrained dynamics 
  /// constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return Squared norm of the residual in the unconstrained dynamics 
  /// constraint.
  ///
  double squaredNormUnconstrainedDynamicsResidual(const double dt) const;

  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Ka,
                            const Eigen::MatrixBase<MatrixType2>& Ku) const;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kaq,
                            const Eigen::MatrixBase<MatrixType2>& Kav,
                            const Eigen::MatrixBase<MatrixType3>& Kuq,
                            const Eigen::MatrixBase<MatrixType4>& Kuv) const;

private:
  Eigen::VectorXd ID_, lu_condensed_;
  Eigen::MatrixXd dID_dq_, dID_dv_, dID_da_, Quu_, Quu_dID_dq_, Quu_dID_dv_, 
                  Quu_dID_da_;
  int dimv_;

  void linearizeInverseDynamics(Robot& robot, const SplitSolution& s);

  void computeInverseDynamicsResidual(Robot& robot, const SplitSolution& s);

};

} // namespace idocp 

#include "idocp/unocp/unconstrained_dynamics.hxx"

#endif // IDOCP_UNCONSTRAINED_DYNAMICS_HPP_ 