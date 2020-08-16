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

  ContactCost(const Robot& robot)
    : CostFunctionComponentBase(),
      f_ref_(Eigen::VectorXd::Zero(robot.max_dimf())),
      f_weight_(Eigen::VectorXd::Zero(robot.max_dimf())) {
  }

  ContactCost()
    : f_ref_(),
      f_weight_() {
  }

  ~ContactCost() {
  }

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
           const double dtau, const SplitSolution& s) const override {
    double l = 0;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (robot.is_contact_active(i)) {
        for (int j=0; j<3; ++j) {
          l += f_weight_.coeff(3*i+j) * (s.f.coeff(3*i+j)-f_ref_.coeff(3*i+j)) 
                                      * (s.f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
        }
      }
    }
    return 0.5 * dtau * l;
  }

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const override {
    return 0;
  }

  template <typename VectorType>
  void lq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lq) const override {
    // do nothing
  }

  template <typename VectorType>
  void lv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lv) const override {
    // do nothing
  }

  template <typename VectorType>
  void la(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& la) const override {
    // do nothing
  }

  template <typename VectorType>
  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          const Eigen::MatrixBase<VectorType>& lf) const override {
    int dimf = 0;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (robot.is_contact_active(i)) {
        for (int j=0; j<3; ++j) {
          const_cast<Eigen::MatrixBase<VectorType>&>(lf).coeffRef(dimf+j) 
              = dtau * f_weight_.coeff(3*i+j) * (s.f.coeff(3*i+j)-f_ref_.coeff(3*i+j));
        }
        dimf += 3;
      }
    }
  }

  template <typename VectorType>
  void lu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<VectorType>& lu) const override {
    // do nothing
  }

  template <typename MatrixType>
  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lqq) const override {
    // do nothing
  }

  template <typename MatrixType>
  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lvv) const override {
    // do nothing
  }

  template <typename MatrixType>
  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& laa) const override {
    // do nothing
  }

  template <typename MatrixType>
  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& lff) const override {
    int dimf = 0;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (robot.is_contact_active(i)) {
        for (int j=0; j<3; ++j) {
          const_cast<Eigen::MatrixBase<MatrixType>&>(lff).coeffRef(dimf+j, dimf+j) 
              = dtau * f_weight_.coeff(3*i+j);
        }
        dimf += 3;
      }
    }
  }

  template <typename MatrixType>
  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           const Eigen::MatrixBase<MatrixType>& luu) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_lqq(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& lqq) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_lvv(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& lvv) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_laa(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& laa) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   const Eigen::MatrixBase<MatrixType>& lff) const override {
    int dimf = 0;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (robot.is_contact_active(i)) {
        for (int j=0; j<3; ++j) {
          const_cast<Eigen::MatrixBase<MatrixType>&>(lff).coeffRef(dimf+j, dimf+j) 
              += dtau * f_weight_.coeff(3*i+j);
        }
        dimf += 3;
      }
    }
  }

  template <typename MatrixType>
  void augment_luu(const Robot& robot, CostFunctionData& data, 
                   const double t, const double dtau, const SplitSolution& s,
                   const Eigen::MatrixBase<MatrixType>& luu) const override {
    // do nothing
  }

  template <typename VectorType>
  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s,
            const Eigen::MatrixBase<VectorType>& phiq) const override {
    // do nothing
  }

  template <typename VectorType>
  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s,
            const Eigen::MatrixBase<VectorType>& phiv) const override {
    // do nothing
  }

  template <typename MatrixType>
  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s,
             const Eigen::MatrixBase<MatrixType>& phiqq) const override {
    // do nothing
  }

  template <typename MatrixType>
  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s,
             const Eigen::MatrixBase<MatrixType>& phivv) const override {
    // do nothing
  }

  template <typename MatrixType>
  void augment_phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phiqq) const override {
    // do nothing
  } 

  template <typename MatrixType>
  void augment_phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s, 
                     const Eigen::MatrixBase<MatrixType>& phivv) const override {
    // do nothing
  } 

private:
  int max_dimf_;
  Eigen::VectorXd f_ref_, f_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_COST_HPP_