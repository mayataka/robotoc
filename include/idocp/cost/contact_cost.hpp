#ifndef IDOCP_CONTACT_COST_HPP_
#define IDOCP_CONTACT_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"


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
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           const Eigen::Ref<const Eigen::VectorXd>& f, 
           const Eigen::Ref<const Eigen::VectorXd>& u) const override;

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
          Eigen::Ref<Eigen::VectorXd> lf) const override;

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
           Eigen::Ref<Eigen::MatrixXd> lff) const override;

  void augment_lff(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& f, 
                   Eigen::Ref<Eigen::MatrixXd> lff) const override;


  // The following functions do nothig, just for dynamic polymorphism.
  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v) const override {
    return 0;
  }

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lq) const override {}

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> lv) const override {}

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
          const Eigen::Ref<const Eigen::VectorXd>& v, 
          const Eigen::Ref<const Eigen::VectorXd>& a, 
          Eigen::Ref<Eigen::VectorXd> la) const override {}

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
          Eigen::Ref<Eigen::VectorXd> lu) const override {}

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lqq) const override {}

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> lvv) const override {}

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
           const Eigen::Ref<const Eigen::VectorXd>& v, 
           const Eigen::Ref<const Eigen::VectorXd>& a, 
           Eigen::Ref<Eigen::MatrixXd> laa) const override {}

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
           Eigen::Ref<Eigen::MatrixXd> luu) const override {}

  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lqq) const override {}

  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lvv) const override {}

  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> laa) const override {}

  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& u, 
                   Eigen::Ref<Eigen::MatrixXd> luu) const override {}

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiq) const override {}

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const Eigen::Ref<const Eigen::VectorXd>& q, 
            const Eigen::Ref<const Eigen::VectorXd>& v, 
            Eigen::Ref<Eigen::VectorXd> phiv) const override {}

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phiqq) const override {}

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const Eigen::Ref<const Eigen::VectorXd>& q, 
             const Eigen::Ref<const Eigen::VectorXd>& v, 
             Eigen::Ref<Eigen::MatrixXd> phivv) const override {}

  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phiqq) const override {}

  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phivv) const override {}

private:
  Eigen::VectorXd f_ref_, f_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_COST_HPP_