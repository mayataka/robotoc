#include "idocp/cost/task_space_6d_cost.hpp"

#include <iostream>


namespace idocp {

TaskSpace6DCost::TaskSpace6DCost(const Robot& robot, const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    SE3_ref_(pinocchio::SE3(Eigen::Matrix3d::Identity(), 
                            Eigen::Vector3d::Zero())),
    SE3_ref_inv_(SE3_ref_.inverse()),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::TaskSpace6DCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    SE3_ref_(pinocchio::SE3(Eigen::Matrix3d::Identity(), 
                            Eigen::Vector3d::Zero())),
    SE3_ref_inv_(SE3_ref_.inverse()),
    q_6d_weight_(Eigen::VectorXd::Zero(6)), 
    qf_6d_weight_(Eigen::VectorXd::Zero(6)) {
}


TaskSpace6DCost::~TaskSpace6DCost() {
}


bool TaskSpace6DCost::useKinematics() const {
  return true;
}


void TaskSpace6DCost::set_q_6d_ref(const pinocchio::SE3& SE3_ref) {
  SE3_ref_ = SE3_ref;
  SE3_ref_inv_ = SE3_ref_.inverse();
}


void TaskSpace6DCost::set_q_6d_ref(const Eigen::Vector3d& position_ref, 
                                   const Eigen::Matrix3d& rotation_mat_ref) {
  SE3_ref_ = pinocchio::SE3(rotation_mat_ref, position_ref);
  SE3_ref_inv_ = SE3_ref_.inverse();
}


void TaskSpace6DCost::set_q_6d_weight(const Eigen::Vector3d& position_weight, 
                                      const Eigen::Vector3d& rotation_weight) {
  q_6d_weight_.template head<3>() = rotation_weight;
  q_6d_weight_.template tail<3>() = position_weight;
}


void TaskSpace6DCost::set_q_6d_weight(const Vector6d& q_6d_weight) {
  q_6d_weight_ = q_6d_weight;
}


void TaskSpace6DCost::set_qf_6d_weight(const Eigen::Vector3d& position_weight, 
                                       const Eigen::Vector3d& rotation_weight) {
  qf_6d_weight_.template head<3>() = rotation_weight;
  qf_6d_weight_.template tail<3>() = position_weight;
}


void TaskSpace6DCost::set_qf_6d_weight(const Vector6d& qf_6d_weight) {
  qf_6d_weight_ = qf_6d_weight;
}


double TaskSpace6DCost::l(Robot& robot, CostFunctionData& data, 
                          const double t, const double dtau, 
                          const SplitSolution& s) const {
  double l = 0;
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  l += (q_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * dtau * l;
}


double TaskSpace6DCost::phi(Robot& robot, CostFunctionData& data, 
                            const double t, const SplitSolution& s) const {
  double phi = 0;
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  phi += (qf_6d_weight_.array()*data.diff_6d.array()*data.diff_6d.array()).sum();
  return 0.5 * phi;
}


void TaskSpace6DCost::lq(Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, 
                         KKTResidual& kkt_residual) const {
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += dtau * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() 
                                       * data.diff_6d;
}


void TaskSpace6DCost::lqq(Robot& robot, CostFunctionData& data, 
                          const double t, const double dtau, 
                          const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_matrix.Qqq().noalias()
      += dtau * data.JJ_6d.transpose() * q_6d_weight_.asDiagonal() * data.JJ_6d;
}


void TaskSpace6DCost::phiq(Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s,
                           KKTResidual& kkt_residual) const {
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  data.diff_6d = pinocchio::log6(data.diff_SE3).toVector();
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_residual.lq().noalias() 
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.diff_6d;
}


void TaskSpace6DCost::phiqq(Robot& robot, CostFunctionData& data, 
                            const double t, const SplitSolution& s,
                            KKTMatrix& kkt_matrix) const {
  data.diff_SE3 = SE3_ref_inv_ * robot.framePlacement(frame_id_);
  pinocchio::Jlog6(data.diff_SE3, data.J_66);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.JJ_6d.noalias() = data.J_66 * data.J_6d;
  kkt_matrix.Qqq().noalias()
      += data.JJ_6d.transpose() * qf_6d_weight_.asDiagonal() * data.JJ_6d;
}

} // namespace idocp