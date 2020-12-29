#ifndef IDOCP_UNCONSTRAINED_DYNAMICS_HPP_
#define IDOCP_UNCONSTRAINED_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


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

  void linearizeUnconstrainedDynamics(Robot& robot, const double dtau, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual);

  template <typename MatrixType1, typename MatrixType2>
  void condenseUnconstrainedDynamics(Robot& robot, const double dtau, 
                                     SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual,
                                     const Eigen::MatrixBase<MatrixType1>& Qaq, 
                                     const Eigen::MatrixBase<MatrixType2>& Qav);
  
  void computeCondensedDirection(const double dtau, 
                                 const SplitKKTMatrix& kkt_matrix, 
                                 const SplitKKTResidual& kkt_residual, 
                                 SplitDirection& d);

  void computeUnconstrainedDynamicsResidual(Robot& robot, 
                                            const SplitSolution& s);

  double l1NormUnconstrainedDynamicsResidual(const double dtau) const;

  double squaredNormUnconstrainedDynamicsResidual(const double dtau) const;

  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Ka,
                            const Eigen::MatrixBase<MatrixType2>& Ku) const;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kaq,
                            const Eigen::MatrixBase<MatrixType2>& Kav,
                            const Eigen::MatrixBase<MatrixType3>& Kuq,
                            const Eigen::MatrixBase<MatrixType4>& Kuv) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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