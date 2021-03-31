#ifndef IDOCP_PERIODIC_COM_REF_HPP_
#define IDOCP_PERIODIC_COM_REF_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/time_varying_com_cost.hpp"


namespace idocp {

class PeriodicCoMRef : public TimeVaryingCoMRefBase {
public:
  PeriodicCoMRef(const Eigen::Vector3d CoM_ref0, const Eigen::Vector3d v_CoM_ref, 
                 const double t0, const double period_swing, 
                 const double period_stance, bool first_step);

  ~PeriodicCoMRef();

  void update_CoM_ref(const double t, Eigen::VectorXd& CoM_ref) const override;

  bool isActive(const double t) const override;

private:
  Eigen::Vector3d CoM_ref0_, v_CoM_ref_;
  double step_length_, t0_, period_swing_, period_stance_, period_;
  bool first_step_;

};

} // namespace idocp


#endif // IDOCP_PERIODIC_COM_REF_HPP_ 