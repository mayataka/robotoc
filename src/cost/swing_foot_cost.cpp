#include "robotoc/cost/swing_foot_cost.hpp"


namespace robotoc {

SwingFootCost::SwingFootCost(const Robot& robot, const int contact_index, 
                             const std::shared_ptr<SwingFootRefBase>& ref) 
  : CostFunctionComponentBase(),
    contact_index_(contact_index), 
    contact_frame_id_(robot.contactFrames()[contact_index]),
    ref_(ref),
    q_3d_weight_(Eigen::Vector3d::Zero()) {
}


SwingFootCost::SwingFootCost()
  : CostFunctionComponentBase(),
    contact_index_(0), 
    contact_frame_id_(0),
    ref_(),
    q_3d_weight_() {
}


SwingFootCost::~SwingFootCost() {
}


void SwingFootCost::set_ref(const std::shared_ptr<SwingFootRefBase>& ref) {
  ref_ = ref;
}


void SwingFootCost::set_q_weight(const Eigen::Vector3d& q_3d_weight) {
  q_3d_weight_ = q_3d_weight;
}


bool SwingFootCost::useKinematics() const {
  return true;
}


double SwingFootCost::evalStageCost(Robot& robot, 
                                    const ContactStatus& contact_status, 
                                    CostFunctionData& data, 
                                    const double t, const double dt, 
                                    const SplitSolution& s) const {
  if (!contact_status.isContactActive(contact_index_)) {
    double l = 0;
    ref_->update_q_3d_ref(contact_status, data.q_3d_ref);
    data.diff_3d = robot.framePosition(contact_frame_id_) - data.q_3d_ref;
    l += (q_3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * dt * l;
  }
  else {
    return 0.0;
  }
}


void SwingFootCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (!contact_status.isContactActive(contact_index_)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(contact_frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(contact_frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += dt * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.diff_3d;
  }
}


void SwingFootCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (!contact_status.isContactActive(contact_index_)) {
    kkt_matrix.Qqq().noalias()
        += dt * data.J_3d.transpose() * q_3d_weight_.asDiagonal() * data.J_3d;
  }
}


double SwingFootCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  return 0.0;
}


void SwingFootCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  // Do nothing
}


void SwingFootCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  // Do nothing
}


double SwingFootCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s) const {
  return 0.0;
}


void SwingFootCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  // Do nothing
}


void SwingFootCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  // Do nothing
}

} // namespace robotoc