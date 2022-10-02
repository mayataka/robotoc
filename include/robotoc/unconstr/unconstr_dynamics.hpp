#ifndef ROBOTOC_UNCONSTR_DYNAMICS_HPP_
#define ROBOTOC_UNCONSTR_DYNAMICS_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"


namespace robotoc {

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
  void evalUnconstrDynamics(Robot& robot, const SplitSolution& s);

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

  template <int p=1>
  double primalFeasibility() const {
    return ID_.template lpNorm<p>();
  }

  ///
  /// @brief Returns squared norm of the KKT residual, that is, the residual in 
  /// the unconstrained dynamics constraint. 
  /// @return Squared norm of the residual in the unconstrained dynamics 
  /// constraint.
  ///
  double KKTError() const {
    return ID_.squaredNorm();
  }

  ///
  /// @brief Returns the lp norm of the constraint violation, that is,
  /// the primal residual in the unconstr dynamics. Default norm is l1-norm.
  /// You can specify l-infty norm by passing Eigen::Infinity as the 
  /// template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the constraint violation.
  ///
  template <int p=1>
  double constraintViolation() const {
    return ID_.template lpNorm<p>();
  }

  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Ka,
                            const Eigen::MatrixBase<MatrixType2>& Ku) const {
    assert(Ka.rows() == dimv_);
    assert(Ka.cols() == 2*dimv_);
    assert(Ku.rows() == dimv_);
    assert(Ku.cols() == 2*dimv_);
    getStateFeedbackGain(
        Ka.leftCols(dimv_), Ka.rightCols(dimv_),
        const_cast<Eigen::MatrixBase<MatrixType2>&>(Ku).leftCols(dimv_),
        const_cast<Eigen::MatrixBase<MatrixType2>&>(Ku).rightCols(dimv_));
  }

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kaq,
                            const Eigen::MatrixBase<MatrixType2>& Kav,
                            const Eigen::MatrixBase<MatrixType3>& Kuq,
                            const Eigen::MatrixBase<MatrixType4>& Kuv) const {
    assert(Kaq.rows() == dimv_);
    assert(Kaq.cols() == dimv_);
    assert(Kav.rows() == dimv_);
    assert(Kav.cols() == dimv_);
    assert(Kuq.rows() == dimv_);
    assert(Kuq.cols() == dimv_);
    assert(Kuv.rows() == dimv_);
    assert(Kuv.cols() == dimv_);
    const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq) = dID_dq_;
    const_cast<Eigen::MatrixBase<MatrixType3>&>(Kuq).noalias() += dID_da_ * Kaq;
    const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv) = dID_dv_;
    const_cast<Eigen::MatrixBase<MatrixType4>&>(Kuv).noalias() += dID_da_ * Kav;
  }

private:
  Eigen::VectorXd ID_, lu_condensed_;
  Eigen::MatrixXd dID_dq_, dID_dv_, dID_da_, Quu_, Quu_dID_dq_, Quu_dID_dv_, 
                  Quu_dID_da_;
  int dimv_;

};

} // namespace robotoc 

#endif // ROBOTOC_UNCONSTR_DYNAMICS_HPP_ 