#ifndef IDOCP_IMPULSE_IMPACT_COST_HPP_
#define IDOCP_IMPULSE_IMPACT_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/impulse_cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

class ImpulseForceCost final : public ImpulseCostFunctionComponentBase {
public:
  ImpulseForceCost(const Robot& robot);

  ImpulseForceCost();

  ~ImpulseForceCost();

  // Use defalut copy constructor.
  ImpulseForceCost(const ImpulseForceCost&) = default;

  // Use defalut copy operator.
  ImpulseForceCost& operator=(const ImpulseForceCost&) = default;

  // Use defalut move constructor.
  ImpulseForceCost(ImpulseForceCost&&) noexcept = default;

  // Use defalut move assign operator.
  ImpulseForceCost& operator=(ImpulseForceCost&&) noexcept = default;

  void set_f_ref(const std::vector<Eigen::Vector3d>& f_ref);

  void set_f_weight(const std::vector<Eigen::Vector3d>& f_weight);

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s) const override;

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override {}

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override {}

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override {}

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override {}

  void ldv(Robot& robot, CostFunctionData& data, const double t, 
           const Eigen::VectorXd& dv, 
           Eigen::VectorXd& ldv) const override {}

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const ImpulseSplitSolution& s, 
          ImpulseKKTResidual& kkt_residual) const override;

  void ldvdv(Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::VectorXd& dv, 
             Eigen::MatrixXd& Qdvdv) const override {}

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const ImpulseSplitSolution& s, 
           ImpulseKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int max_point_contacts_, max_dimf_;
  std::vector<Eigen::Vector3d> f_ref_, f_weight_;

};

} // namespace idocp


#endif // IDOCP_IMPULSE_IMPACT_COST_HPP_ 