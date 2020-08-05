#ifndef IDOCP_JOINT_SPACE_COST_HPP_
#define IDOCP_JOINT_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class JointSpaceCost {
public:
  JointSpaceCost(const Robot& robot, const Eigen::VectorXd& q_weight,  
                 const Eigen::VectorXd& v_weight, 
                 const Eigen::VectorXd& a_weight,  
                 const Eigen::VectorXd& u_weight,
                 const Eigen::VectorXd& qf_weight,  
                 const Eigen::VectorXd& vf_weight);

  JointSpaceCost(const Robot& robot, const Eigen::VectorXd& q_ref,  
                 const Eigen::VectorXd& v_ref, const Eigen::VectorXd& a_ref,  
                 const Eigen::VectorXd& u_ref, const Eigen::VectorXd& q_weight,  
                 const Eigen::VectorXd& v_weight, 
                 const Eigen::VectorXd& a_weight,  
                 const Eigen::VectorXd& u_weight,
                 const Eigen::VectorXd& qf_weight,  
                 const Eigen::VectorXd& vf_weight);

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

  double l(const double dtau, const Eigen::VectorXd& q, 
           const Eigen::VectorXd& v, const Eigen::VectorXd& a,
           const Eigen::VectorXd& u) const;

  void lq(const Robot& robot, const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& lq);

  void lv(const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& lv) const;

  void la(const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& la) const;

  void lu(const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& lu) const;

  void lqq(const Robot& robot, const double dtau, Eigen::MatrixXd& lqq) const;

  void lvv(const double dtau, Eigen::MatrixXd& lvv) const;

  void laa(const double dtau, Eigen::MatrixXd& laa) const;

  void luu(const double dtau, Eigen::MatrixXd& luu) const;

  void augment_lqq(const Robot& robot, const double dtau, 
                  Eigen::MatrixXd& lqq) const;

  void augment_lvv(const double dtau, Eigen::MatrixXd& lvv) const;

  void augment_laa(const double dtau, Eigen::MatrixXd& laa) const;

  void augment_luu(const double dtau, Eigen::MatrixXd& luu) const;

  double phi(const Eigen::VectorXd& q, const Eigen::VectorXd& v) const;

  void phiq(const Robot& robot, const Eigen::VectorXd& q, Eigen::VectorXd& phiq);

  void phiv(const Eigen::VectorXd& v, Eigen::VectorXd& phiv) const;

  void phiqq(const Robot& robot, Eigen::MatrixXd& phiqq) const;

  void phivv(Eigen::MatrixXd& phivv) const;

private:
  bool has_floating_base_;
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_, 
                  lq_configuration_;
  Eigen::MatrixXd lqq_configuration_, phiqq_configuration_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_HPP_