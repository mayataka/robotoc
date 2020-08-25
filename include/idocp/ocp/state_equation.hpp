#ifndef IDOCP_STATE_EQUATION_HPP_
#define IDOCP_STATE_EQUATION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class StateEquation {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  StateEquation(const Robot& robot);

  StateEquation();

  ~StateEquation();

  StateEquation(const StateEquation&) = default;

  StateEquation& operator=(const StateEquation&) = default;
 
  StateEquation(StateEquation&&) noexcept = default;

  StateEquation& operator=(StateEquation&&) noexcept = default;

  void linearizeForwardEuler(Robot& robot, const double dtau, 
                             const SplitSolution& s, 
                             const SplitSolution& s_next, 
                             KKTMatrix& kkt_matrix, 
                             KKTResidual& kkt_residual) const;

  template <typename ConfigVectorType1, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3, 
            typename ConfigVectorType2>
  void linearizeBackwardEuler(
      Robot& robot, const double dtau, 
      const Eigen::MatrixBase<ConfigVectorType1>& q_prev, 
      const Eigen::MatrixBase<TangentVectorType1>& v_prev, 
      const SplitSolution& s, 
      const Eigen::MatrixBase<TangentVectorType2>& lmd_next, 
      const Eigen::MatrixBase<TangentVectorType3>& gmm_next, 
      const Eigen::MatrixBase<ConfigVectorType2>& q_next, 
      KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  template <typename ConfigVectorType, typename TangentVectorType>
  void linearizeBackwardEulerTerminal(
      Robot& robot, const double dtau, 
      const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const Eigen::MatrixBase<TangentVectorType>& v_prev, 
      const SplitSolution& s, 
      KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) const;

  double violationL1Norm(const KKTResidual& kkt_residual) const;

  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3>
  double computeForwardEulerViolationL1Norm(
      const Robot& robot, const double step_size, const double dtau, 
      const SplitSolution& s, const Eigen::MatrixBase<ConfigVectorType>& q_next, 
      const Eigen::MatrixBase<TangentVectorType1>& v_next, 
      const Eigen::MatrixBase<TangentVectorType2>& dq_next, 
      const Eigen::MatrixBase<TangentVectorType3>& dv_next, 
      KKTResidual& kkt_residual) const;

  template <typename ConfigVectorType, typename TangentVectorType>
  double computeBackwardEulerViolationL1Norm(
      const Robot& robot, const double dtau, 
      const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const Eigen::MatrixBase<TangentVectorType>& v_prev, 
      const SplitSolution& s, KKTResidual& kkt_residual) const;

  template <typename ConfigVectorType, typename TangentVectorType1, 
            typename TangentVectorType2, typename TangentVectorType3>
  double computeBackwardEulerViolationL1Norm(
      const Robot& robot, const double step_size, const double dtau, 
      const Eigen::MatrixBase<ConfigVectorType>& q_prev, 
      const Eigen::MatrixBase<TangentVectorType1>& v_prev, 
      const Eigen::MatrixBase<TangentVectorType2>& dq_prev, 
      const Eigen::MatrixBase<TangentVectorType3>& dv_prev, 
      const SplitSolution& s, KKTResidual& kkt_residual) const;

private:
  Eigen::MatrixXd dsubtract_dq_;

};

} // namespace idocp 

#include "idocp/ocp/state_equation.hxx"

#endif // IDOCP_STATE_EQUATION_HPP_