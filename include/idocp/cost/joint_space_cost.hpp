#ifndef IDOCP_JOINT_SPACE_COST_HPP_
#define IDOCP_JOINT_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class JointSpaceCost final : public CostFunctionComponentBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  JointSpaceCost(const Robot& robot);

  JointSpaceCost();

  ~JointSpaceCost();

  // Use defalut copy constructor.
  JointSpaceCost(const JointSpaceCost&) = default;

  // Use defalut copy operator.
  JointSpaceCost& operator=(const JointSpaceCost&) = default;

  // Use defalut move constructor.
  JointSpaceCost(JointSpaceCost&&) noexcept = default;

  // Use defalut move assign operator.
  JointSpaceCost& operator=(JointSpaceCost&&) noexcept = default;

  void set_q_ref(const Eigen::VectorXd& q_ref);

  void set_v_ref(const Eigen::VectorXd& v_ref);

  void set_a_ref(const Eigen::VectorXd& a_ref);
 
  void set_u_ref(const Eigen::VectorXd& u_ref);

  void set_q_weight(const Eigen::VectorXd& q_weight);

  void set_v_weight(const Eigen::VectorXd& v_weight);

  void set_a_weight(const Eigen::VectorXd& a_weight);

  void set_u_weight(const Eigen::VectorXd& u_weight);

  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           const Eigen::Ref<const Eigen::VectorXd>& f, 
           const Eigen::Ref<const Eigen::VectorXd>& u) const override;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v) const override;

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lq) const override;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lv) const override;

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> la) const override;

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
          Eigen::Ref<Eigen::VectorXd> lu) const override;

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lqq) const override;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lvv) const override;

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> laa) const override;

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
           Eigen::Ref<Eigen::MatrixXd> luu) const override;

  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lqq) const override;

  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lvv) const override;

  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> laa) const override;

  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& u, 
                   Eigen::Ref<Eigen::MatrixXd> luu) const override;

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiq) const override;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiv) const override;

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phiqq) const override;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phivv) const override;

  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phiqq) const override;

  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phivv) const override;


  // The following functions do nothig, just for dynamic polymorphism.
  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
          Eigen::Ref<Eigen::VectorXd> lf) const override {}

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
           Eigen::Ref<Eigen::MatrixXd> lff) const override {}

  void augment_lff(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& f, 
                   Eigen::Ref<Eigen::MatrixXd> lff) const override {}
private:
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_HPP_