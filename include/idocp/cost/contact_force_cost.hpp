#ifndef IDOCP_CONTACT_FORCE_COST_HPP_
#define IDOCP_CONTACT_FORCE_COST_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class ContactForceCost final : public CostFunctionComponentBase {
public:
  ContactForceCost(const Robot& robot);

  ContactForceCost();

  ~ContactForceCost();

  // Use defalut copy constructor.
  ContactForceCost(const ContactForceCost&) = default;

  // Use defalut copy operator.
  ContactForceCost& operator=(const ContactForceCost&) = default;

  // Use defalut move constructor.
  ContactForceCost(ContactForceCost&&) noexcept = default;

  // Use defalut move assign operator.
  ContactForceCost& operator=(ContactForceCost&&) noexcept = default;

  bool useKinematics() const override;

  void set_f_ref(const std::vector<Eigen::Vector3d>& f_ref);

  void set_f_ref(const Robot& robot);

  void set_f_weight(const std::vector<Eigen::Vector3d>& f_weight);

  void set_fi_ref(const std::vector<Eigen::Vector3d>& fi_ref);

  void set_fi_weight(const std::vector<Eigen::Vector3d>& fi_weight);

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override;

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override {}

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const override {}

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override;


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int max_point_contacts_, max_dimf_;
  std::vector<Eigen::Vector3d> f_ref_, f_weight_, fi_ref_, fi_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_FORCE_COST_HPP_ 