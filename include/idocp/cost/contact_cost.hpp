#ifndef IDOCP_CONTACT_COST_HPP_
#define IDOCP_CONTACT_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class ContactCost final : public CostFunctionComponentBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ContactCost(const Robot& robot);

  ContactCost();

  ~ContactCost();

  // Use defalut copy constructor.
  ContactCost(const ContactCost&) = default;

  // Use defalut copy operator.
  ContactCost& operator=(const ContactCost&) = default;

  // Use defalut move constructor.
  ContactCost(ContactCost&&) noexcept = default;

  // Use defalut move assign operator.
  ContactCost& operator=(ContactCost&&) noexcept = default;

  void set_f_ref(const Eigen::VectorXd& f_ref);

  void set_f_weight(const Eigen::VectorXd& f_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const override;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override; 

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override {}

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s,
          KKTResidual& kkt_residual) const override {}

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const override;

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override {}

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const override;

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override {}

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const override {}

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override {}

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const override {}

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const override {}

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const override {}

private:
  int max_dimf_;
  Eigen::VectorXd f_ref_, f_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_COST_HPP_