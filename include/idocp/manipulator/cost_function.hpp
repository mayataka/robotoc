#ifndef IDOCP_MANIPULATOR_COST_FUNCTION_HPP_
#define IDOCP_MANIPULATOR_COST_FUNCTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"


namespace idocp {
namespace manipulator {

class CostFunction final : public CostFunctionInterface {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CostFunction(const Robot& robot);

  CostFunction();

  ~CostFunction();

  CostFunction(const CostFunction&) = default;

  CostFunction& operator=(const CostFunction&) = default;

  CostFunction(CostFunction&&) noexcept = default;

  CostFunction& operator=(CostFunction&&) noexcept = default;

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_a_ref(const Eigen::VectorXd& a_ref);

  void set_u_ref(const Eigen::VectorXd& u_ref);

  void set_f_ref(const Eigen::VectorXd& f_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_u_weight(const Eigen::VectorXd& u_weight);

  void set_f_weight(const Eigen::VectorXd& f_weight);

  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           const Eigen::VectorXd& u, const Eigen::VectorXd& f) const override;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::VectorXd& q, 
             const Eigen::VectorXd& v) const override;

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& lq) const override;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& lv) const override;

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& la) const override;

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& la) const override;

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& la) const override;

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& lqq) const override;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& lvv) const override;

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& laa) const override;

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& luu) const override;

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& f, 
           Eigen::MatrixXd& lff) const override;

  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lqq) const override;

  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lvv) const override;

  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& laa) const override;

  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& luu) const override;

  void augment_lff(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::VectorXd& f, 
                   Eigen::MatrixXd& lff) const override;

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiq) const override;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiv) const override;

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             Eigen::MatrixXd& phiqq) const override;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             Eigen::MatrixXd& phivv) const override;

private:
  JointSpaceCost joint_space_cost_;  
  ContactCost contact_cost_;

};

} // namespace manipulator
} // namespace idocp


#endif // IDOCP_MANIPULATOR_COST_FUNCTION_HPP_