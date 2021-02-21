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

class UnconstrainedDynamics {
public:
  UnconstrainedDynamics(const Robot& robot);

  UnconstrainedDynamics();

  ~UnconstrainedDynamics();

  UnconstrainedDynamics(const UnconstrainedDynamics&) = default;

  UnconstrainedDynamics& operator=(const UnconstrainedDynamics&) = default;
 
  UnconstrainedDynamics(UnconstrainedDynamics&&) noexcept = default;

  UnconstrainedDynamics& operator=(UnconstrainedDynamics&&) noexcept = default;

  void linearizeUnconstrainedDynamics(Robot& robot, const double dt, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual);

  void condenseUnconstrainedDynamics(const SplitKKTMatrix& kkt_matrix, 
                                     const SplitKKTResidual& kkt_residual,
                                     SplitUnKKTMatrix& unkkt_matrix, 
                                     SplitUnKKTResidual& unkkt_residual);

  void computeCondensedDirection(const double dt, 
                                 const SplitKKTMatrix& kkt_matrix, 
                                 const SplitKKTResidual& kkt_residual, 
                                 SplitDirection& d);

  void computeUnconstrainedDynamicsResidual(Robot& robot, 
                                            const SplitSolution& s);

  double l1NormUnconstrainedDynamicsResidual(const double dt) const;

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