#include "idocp/cost/trotting_foot_step_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TrottingFootStepCost::TrottingFootStepCost(const Robot& robot, 
                                           const int frame_id)
  : CostFunctionComponentBase(),
    frame_id_(frame_id),
    t_start_(0),
    t_period_(0),
    step_length_(0),
    step_height_(0),
    q_3d_ref_init_(Eigen::Vector3d::Zero()),
    q_3d_weight_(Eigen::Vector3d::Zero()),
    qf_3d_weight_(Eigen::Vector3d::Zero()) {
}


TrottingFootStepCost::TrottingFootStepCost()
  : CostFunctionComponentBase(),
    frame_id_(0),
    t_start_(0),
    t_period_(0),
    step_length_(0),
    step_height_(0),
    q_3d_ref_init_(),
    q_3d_weight_(),
    qf_3d_weight_() {
}


TrottingFootStepCost::~TrottingFootStepCost() {
}


bool TrottingFootStepCost::useKinematics() const {
  return true;
}


void TrottingFootStepCost::set_q_3d_ref(const Eigen::Vector3d& q_3d_ref_init,
                                        const double step_length, 
                                        const double step_height) {
  try {
    if (step_length <= 0) {
      throw std::invalid_argument(
          "invalid argument: step_length must be positive!");
    }
    if (step_height < step_height) {
      throw std::runtime_error(
          "invalid argument: step_height must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  q_3d_ref_init_ = q_3d_ref_init;
  step_length_ = step_length;
  step_height_ = step_height;
}


void TrottingFootStepCost::set_period(const double t_start, 
                                      const double t_period) {
  try {
    if (t_period <= 0) {
      throw std::invalid_argument(
          "invalid argument: t_period must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  t_start_ = t_start;
  t_period_ = t_period;
}


void TrottingFootStepCost::set_q_3d_weight(const Eigen::Vector3d& q_3d_weight) {
  q_3d_weight_ = q_3d_weight;
}


void TrottingFootStepCost::set_qf_3d_weight(
    const Eigen::Vector3d& qf_3d_weight) {
  qf_3d_weight_ = qf_3d_weight;
}


double TrottingFootStepCost::l(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s) const {
  double l = 0;
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref(t);
  l += (q_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * dtau * l;
}


double TrottingFootStepCost::phi(Robot& robot, CostFunctionData& data, 
                                 const double t, const SplitSolution& s) const {
  double phi = 0;
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref(t);
  phi += (qf_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
  return 0.5 * phi;
}


void TrottingFootStepCost::lq(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s, 
                              SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref(t);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += dtau * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.diff_3d;
}


void TrottingFootStepCost::lqq(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const {
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_matrix.Qqq().noalias()
      += dtau * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.J_3d;
}


void TrottingFootStepCost::phiq(Robot& robot, CostFunctionData& data, 
                                const double t, const SplitSolution& s,
                                SplitKKTResidual& kkt_residual) const {
  data.diff_3d = robot.framePosition(frame_id_) - q_3d_ref(t);
  robot.getFrameJacobian(frame_id_, data.J_6d);
  data.J_3d.noalias() 
      = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
  kkt_residual.lq().noalias() 
      += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.diff_3d;
}


void TrottingFootStepCost::phiqq(Robot& robot, CostFunctionData& data, 
                                 const double t, const SplitSolution& s,
                                 SplitKKTMatrix& kkt_matrix) const {
    robot.getFrameJacobian(frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(frame_id_) * data.J_6d.template topRows<3>();
    kkt_matrix.Qqq().noalias()
        += data.J_3d.transpose() * qf_3d_weight_.asDiagonal() * data.J_3d;
}

} // namespace idocp