#include "idocp/cost/time_varying_task_space_3d_cost.hpp"

#include <iostream>


namespace idocp {

TimeVaryingTaskSpace3DCost::TimeVaryingTaskSpace3DCost(const Robot& robot, 
                                                       const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    t0_(0),
    q0_(Eigen::Vector3d::Zero()),
    v0_(Eigen::Vector3d::Zero()),
    q_3d_weight_(Eigen::Vector3d::Zero()),
    qf_3d_weight_(Eigen::Vector3d::Zero()) {
}


TimeVaryingTaskSpace3DCost::TimeVaryingTaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(),
    t0_(0),
    q0_(),
    v0_(),
    q_3d_weight_(),
    qf_3d_weight_() {
}


TimeVaryingTaskSpace3DCost::~TimeVaryingTaskSpace3DCost() {
}


bool TimeVaryingTaskSpace3DCost::useKinematics() const {
  return true;
}


void TimeVaryingTaskSpace3DCost::set_ref(const double t0, 
                                         const Eigen::Vector3d& q0, 
                                         const Eigen::Vector3d& v0) {
  t0_ = t0;
  q0_ = q0;
  v0_ = v0;
}


void TimeVaryingTaskSpace3DCost::set_q_3d_weight(
    const Eigen::Vector3d& q_3d_weight) {
  q_3d_weight_ = q_3d_weight;
}


void TimeVaryingTaskSpace3DCost::set_qf_3d_weight(
    const Eigen::Vector3d& qf_3d_weight) {
  qf_3d_weight_ = qf_3d_weight;
}


double TimeVaryingTaskSpace3DCost::l(Robot& robot, CostFunctionData& data, 
                                     const double t, const double dtau, 
                                     const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - (q0_+(t-t0_)*v0_);
  l += (q_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * dtau * l;
}


double TimeVaryingTaskSpace3DCost::phi(Robot& robot, CostFunctionData& data, 
                                       const double t, 
                                       const SplitSolution& s) const {
  double phi = 0;
  data.diff_3d = robot.framePosition(frame_id_) - (q0_+(t-t0_)*v0_);
  phi += (qf_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * phi;
}


void TimeVaryingTaskSpace3DCost::lq(Robot& robot, CostFunctionData& data, 
                                    const double t, const double dtau, 
                                    const SplitSolution& s, 
                                    SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.framePosition(frame_id_) - (q0_+(t-t0_)*v0_);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += dtau * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.diff_3d;
}


void TimeVaryingTaskSpace3DCost::lqq(Robot& robot, CostFunctionData& data, 
                                     const double t, const double dtau, 
                                     const SplitSolution& s, 
                                     SplitKKTMatrix& kkt_matrix) const {
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_matrix.Qqq().noalias()
      += dtau * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.J_3d;
}


void TimeVaryingTaskSpace3DCost::phiq(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s,
                                      SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.framePosition(frame_id_) - (q0_+(t-t0_)*v0_);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.diff_3d;
}


void TimeVaryingTaskSpace3DCost::phiqq(Robot& robot, CostFunctionData& data, 
                                       const double t, const SplitSolution& s,
                                       SplitKKTMatrix& kkt_matrix) const {
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.J_3d;
}

} // namespace idocp