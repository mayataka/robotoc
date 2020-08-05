#ifndef IDOCP_CONTACT_COST_HPP_
#define IDOCP_CONTACT_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class ContactCost {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  ContactCost(const Robot& robot, const Eigen::VectorXd& f_weight);

  ContactCost(const Robot& robot, const Eigen::VectorXd& f_ref, 
              const Eigen::VectorXd& f_weight);

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

  double l(const Robot& robot, const double dtau, 
           const Eigen::VectorXd& f) const;

  void lf(const Robot& robot, const double dtau, const Eigen::VectorXd& f, 
          Eigen::VectorXd& lf) const;

  void lff(const Robot& robot, const double dtau, Eigen::MatrixXd& lff) const;

  void augment_lff(const Robot& robot, const double dtau, 
                   Eigen::MatrixXd& lff) const;

private:
  int max_point_contacts_, max_dimf_;
  Eigen::VectorXd f_ref_, f_weight_;

};

} // namespace idocp


#endif // IDOCP_CONTACT_COST_HPP_