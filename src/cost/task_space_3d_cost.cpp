#include "idocp/cost/task_space_3d_cost.hpp"

#include <iostream>


namespace idocp {

TaskSpace3DCost::TaskSpace3DCost(const Robot& robot, 
                                 const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    q_ref_(Eigen::Vector3d::Zero(robot.dimq())),
    q_weight_(Eigen::Vector3d::Zero(robot.dimv())),
    qf_weight_(Eigen::Vector3d::Zero(robot.dimv())) {
}


TaskSpace3DCost::TaskSpace3DCost()
  : CostFunctionComponentBase(),
    frame_id_(),
    q_ref_(),
    q_weight_(),
    qf_weight_() {
}


TaskSpace3DCost::~TaskSpace3DCost() {
}


void TaskSpace3DCost::set_q_ref(const Eigen::Vector3d& q_ref) {
  q_ref_ = q_ref;
}


void TaskSpace3DCost::set_q_weight(const Eigen::Vector3d& q_weight) {
  q_weight_ = q_weight;
}


void TaskSpace3DCost::set_qf_weight(const Eigen::Vector3d& qf_weight) {
  qf_weight_ = qf_weight;
}


double TaskSpace3DCost::l(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s) const {
  double l = 0;
  data.qdiff_3d = robot.framePosition(frame_id_) - q_ref_;
  l += (q_weight_.array()*data.qdiff_3d.array()*data.qdiff_3d.array.array()).sum();
  return 0.5 * dtau * l;
}


double TaskSpace3DCost::phi(const Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s) const {
  double phi = 0;
  data.qdiff_3d = robot.framePosition(frame_id_) - q_ref_;
  phi += (qf_weight_.array()*data.qdiff_3d.array()*data.qdiff_3d.array.array()).sum();
  return 0.5 * phi;
}


void TaskSpace3DCost::lq(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  data.qdiff_3d = robot.framePosition(frame_id_) - q_ref_;
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += dtau * data.J_3d.transpose() * q_weight_.asDiagonal() * data.qdiff_3d;
}


void TaskSpace3DCost::lqq(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_matrix.Qqq().noalias()
      += dtau * data.J_3d.transpose() * q_weight_.asDiagonal() * data.J_3d;
}


void TaskSpace3DCost::phiq(const Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s,
                          KKTResidual& kkt_residual) const {
  data.qdiff_3d = robot.framePosition(frame_id_) - q_ref_;
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.qdiff_3d;
}


void TaskSpace3DCost::phiqq(const Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s,
                           KKTMatrix& kkt_matrix) const {
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * qf_weight_.asDiagonal() * data.J_3d;
}

} // namespace idocp