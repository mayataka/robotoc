#ifndef IDOCP_QUADRUPED_COST_FUNCTION_HPP_
#define IDOCP_QUADRUPED_COST_FUNCTION_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "cost/cost_function_interface.hpp"
#include "cost/joint_space_cost.hpp"
#include "cost/contact_cost.hpp"


namespace idocp {
namespace quadruped {

class CostFunction final : public CostFunctionInterface {
public:
  CostFunction(const Robot& robot);

  ~CostFunction();

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

  void setConfigurationJacobian(const Robot& robot, 
                                const Eigen::VectorXd& q) override;

  double l(const double t, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           const Eigen::VectorXd& u) override; 

  double phi(const double t, const Eigen::VectorXd& q, 
             const Eigen::VectorXd& v) override;

  void lq(const double t, const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& lq) override;

  void lv(const double t, const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& lv) override;

  void la(const double t, const double dtau, const Eigen::VectorXd& q, 
          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
          Eigen::VectorXd& la) override;

  void lu(const double t, const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) override;

  void lf(const double t, const double dtau, const Eigen::VectorXd& f, 
          Eigen::VectorXd& lf) override;

  void lqq(const double t, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& lqq) override;

  void lvv(const double t, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& lvv) override;

  void laa(const double t, const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
           Eigen::MatrixXd& laa) override;

  void luu(const double t, const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& luu) override;

  void lff(const double t, const double dtau, const Eigen::VectorXd& f, 
           Eigen::MatrixXd& lff) override;

  void augment_lqq(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lqq) override;

  void augment_lvv(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lvv) override;

  void augment_laa(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& laa) override;

  void augment_luu(const double t, const double dtau, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& luu) override;

  void augment_lff(const double t, const double dtau, const Eigen::VectorXd& f, 
                   Eigen::MatrixXd& lff) override;

  void phiq(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiq) override;

  void phiv(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiv) override;

  void phiqq(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             Eigen::MatrixXd& phiqq) override;

  void phivv(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             Eigen::MatrixXd& phivv) override;

private:
  JointSpaceCost joint_space_cost_;  
  ContactCost contact_cost_;

};

} // namespace quadruped
} // namespace idocp

#endif // IDOCP_QUADRUPED_COST_FUNCTION_HPP_