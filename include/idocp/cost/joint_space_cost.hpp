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
           const double dtau, const SplitSolution& s) const override {
    double l = 0;
    if (robot.has_floating_base()) {
      robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
      l += (q_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
    }
    else {
      l += (q_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
    }
    l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
    l += (a_weight_.array()*(s.a-a_ref_).array()*(s.a-a_ref_).array()).sum();
    l += (u_weight_.array()*(s.u-u_ref_).array()*(s.u-u_ref_).array()).sum();
    return 0.5 * dtau * l;
  }

  // The following functions do nothig, just for dynamic polymorphism.
  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override {
    double phi = 0;
    if (robot.has_floating_base()) {
      robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
      phi += (qf_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
    }
    else {
      phi += (qf_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
    }
    phi += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
    return 0.5 * phi;
  }

  template <typename VectorType>
  void lq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lq) const override {
    // do nothing
    if (robot.has_floating_base()) {
      robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      const_cast<Eigen::MatrixBase<VectorType>&>(lq)
          = dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.q_diff;
    }
    else {
      (const_cast<Eigen::MatrixBase<VectorType>&>(lq)).array()
          = dtau * q_weight_.array() * (s.q.array()-q_ref_.array());
    }
  }

  template <typename VectorType>
  void lv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lv) const override {
    (const_cast<Eigen::MatrixBase<VectorType>&>(lv)).array()
        = dtau * v_weight_.array() * (s.v.array()-v_ref_.array());
  }

  template <typename VectorType>
  void la(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& la) const override {
    (const_cast<Eigen::MatrixBase<VectorType>&>(la)).array()
        = dtau * a_weight_.array() * (s.a.array()-a_ref_.array());
  }

  template <typename VectorType>
  void lf(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lq) const override {
    // do nothing
  }

  template <typename VectorType>
  void lu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lu) const override {
    (const_cast<Eigen::MatrixBase<VectorType>&>(lu)).array()
        = dtau * u_weight_.array() * (s.u.array()-u_ref_.array());
  }

  template <typename MatrixType>
  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lqq) const override {
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      const_cast<Eigen::MatrixBase<MatrixType>&>(lqq)
          = dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.Jq_diff;
    }
    else {
      const_cast<Eigen::MatrixBase<MatrixType>&>(lqq)
          = dtau * q_weight_.asDiagonal();
    }
  }

  template <typename MatrixType>
  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lvv) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(lvv)
        = dtau * v_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& laa) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(laa)
        = dtau * a_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lff) const override {
    // do nothing
  }

  template <typename MatrixType>
  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& luu) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(luu)
        = dtau * u_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& lqq) const override {
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      const_cast<Eigen::MatrixBase<MatrixType>&>(lqq).noalias()
          += dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.Jq_diff;
    }
    else {
      const_cast<Eigen::MatrixBase<MatrixType>&>(lqq).noalias()
          += dtau * q_weight_.asDiagonal();
    }
  }

  template <typename MatrixType>
  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& lvv) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(lvv).noalias()
        += dtau * v_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& laa) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(laa).noalias()
        += dtau * a_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void augment_lff(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& lff) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& luu) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(luu).noalias()
        += dtau * u_weight_.asDiagonal();
  }

  template <typename VectorType>
  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s,
            const Eigen::MatrixBase<VectorType>& phiq) const override {
    if (robot.has_floating_base()) {
      robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      const_cast<Eigen::MatrixBase<VectorType>&>(phiq)
          = data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.q_diff;
    }
    else {
      (const_cast<Eigen::MatrixBase<VectorType>&>(phiq)).array()
          = qf_weight_.array() * (s.q.array()-q_ref_.array());
    }
  }

  template <typename VectorType>
  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s,
            const Eigen::MatrixBase<VectorType>& phiv) const override {
    (const_cast<Eigen::MatrixBase<VectorType>&>(phiv)).array()
        = vf_weight_.array() * (s.v.array()-v_ref_.array());
  }

  template <typename MatrixType>
  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s,
             const Eigen::MatrixBase<MatrixType>& phiqq) const override {
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq)
          = data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.Jq_diff;
    }
    else {
      const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq)
          = dtau * qf_weight_.asDiagonal();
    }
  }

  template <typename MatrixType>
  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s,
             const Eigen::MatrixBase<MatrixType>& phivv) const override {
    const_cast<Eigen::MatrixBase<MatrixType>&>(phivv)
        = dtau * vf_weight_.asDiagonal();
  }

  template <typename MatrixType>
  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phiqq) const override {
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      (const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq)).noalias()
          += data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.Jq_diff;
    }
    else {
      (const_cast<Eigen::MatrixBase<MatrixType>&>(phiqq)).noalias()
          += dtau * qf_weight_.asDiagonal();
    }
  } 

  template <typename MatrixType>
  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phivv) const override {
    (const_cast<Eigen::MatrixBase<MatrixType>&>(phivv)).noalias()
        += dtau * vf_weight_.asDiagonal();
  } 


private:
  int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_, qf_weight_, vf_weight_;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_HPP_