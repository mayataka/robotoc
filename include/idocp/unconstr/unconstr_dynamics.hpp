#ifndef IDOCP_UNCONSTR_DYNAMICS_HPP_
#define IDOCP_UNCONSTR_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class UnconstrDynamics
/// @brief Inverse dynamics constraint without constraints (floating base or 
/// contacts).
///
class UnconstrDynamics {
public:
  ///
  /// @brief Constructs the inverse dynamics.
  /// @param[in] robot Robot model. 
  ///
  UnconstrDynamics(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrDynamics();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrDynamics();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrDynamics(const UnconstrDynamics&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrDynamics& operator=(const UnconstrDynamics&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrDynamics(UnconstrDynamics&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrDynamics& operator=(UnconstrDynamics&&) noexcept = default;

  ///
  /// @brief Computes the residual in the unconstrained dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this time stage.
  ///
  void computeUnconstrDynamicsResidual(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Linearizes the unconstrained dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeUnconstrDynamics(Robot& robot, const double dt, 
                                 const SplitSolution& s, 
                                 SplitKKTResidual& kkt_residual);

  ///
  /// @brief Condenses the unconstrained dynamics constraint. 
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  ///
  void condenseUnconstrDynamics(SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual);

  ///
  /// @brief Expands the primal variables, i.e., computes the Newton direction 
  /// of the condensed primal variable (control input torques) of this stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void expandPrimal(SplitDirection& d) const;

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variable (Lagrange multiplier) of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  static void expandDual(const double dt, const SplitKKTMatrix& kkt_matrix, 
                         const SplitKKTResidual& kkt_residual, 
                         SplitDirection& d);

  ///
  /// @brief Returns l1-norm of the residual in the unconstrained dynamics 
  /// constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return l1-norm of the residual in the unconstrained dynamics constraint.
  ///
  double l1NormUnconstrDynamicsResidual(const double dt) const;

  ///
  /// @brief Returns squared norm of the residual in the unconstrained dynamics 
  /// constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return Squared norm of the residual in the unconstrained dynamics 
  /// constraint.
  ///
  double squaredNormUnconstrDynamicsResidual(const double dt) const;

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

};

} // namespace idocp 

#include "idocp/unconstr/unconstr_dynamics.hxx"

#endif // IDOCP_UNCONSTR_DYNAMICS_HPP_ 