#ifndef IDOCP_PERIODIC_FOOT_TRACK_REF_HPP_
#define IDOCP_PERIODIC_FOOT_TRACK_REF_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/time_varying_task_space_3d_cost.hpp"


namespace idocp {

class PeriodicFootTrackRef : public TimeVaryingTaskSpace3DRefBase {
public:
  PeriodicFootTrackRef(const Eigen::Vector3d p0, const double step_length, 
                       const double step_height, const double t0, 
                       const double period_swing, const double period_stance, 
                       bool first_step);

  ~PeriodicFootTrackRef();

  void update_q_3d_ref(const double t, Eigen::VectorXd& q_3d_ref) const override;

  bool isActive(const double t) const override;

private:
  Eigen::Vector3d p0_;
  double step_length_, step_height_, t0_, period_swing_, period_stance_,
         period_;
  bool first_step_;

};

} // namespace idocp


#endif // IDOCP_PERIODIC_FOOT_TRACK_REF_HPP_ 