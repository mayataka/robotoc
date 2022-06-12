#include "robotoc/cost/swing_foot_cost.hpp"


namespace robotoc {

SwingFootCost::SwingFootCost(const Robot& robot, const int contact_index, 
                             const std::shared_ptr<SwingFootRefBase>& x3d_ref) 
  : CostFunctionComponentBase(),
    contact_index_(contact_index), 
    contact_frame_id_(robot.contactFrames()[contact_index]),
    x3d_ref_(x3d_ref),
    x3d_weight_(Eigen::Vector3d::Zero()) {
}


SwingFootCost::SwingFootCost(const Robot& robot, 
                             const std::string& contact_frame_name,
                             const std::shared_ptr<SwingFootRefBase>& x3d_ref) 
  : SwingFootCost(robot, 
                  robot.createContactStatus().findContactIndex(contact_frame_name),
                  x3d_ref) {
}


SwingFootCost::SwingFootCost()
  : CostFunctionComponentBase(),
    contact_index_(0), 
    contact_frame_id_(0),
    x3d_ref_(),
    x3d_weight_() {
}


SwingFootCost::~SwingFootCost() {
}


void SwingFootCost::set_x3d_ref(
    const std::shared_ptr<SwingFootRefBase>& x3d_ref) {
  x3d_ref_ = x3d_ref;
}


void SwingFootCost::set_x3d_weight(const Eigen::Vector3d& x3d_weight) {
  x3d_weight_ = x3d_weight;
}


bool SwingFootCost::useKinematics() const {
  return true;
}


double SwingFootCost::evalStageCost(Robot& robot, 
                                    const ContactStatus& contact_status, 
                                    CostFunctionData& data, 
                                    const GridInfo& grid_info, 
                                    const SplitSolution& s) const {
  if (!contact_status.isContactActive(contact_index_)) {
    double l = 0;
    x3d_ref_->update_x3d_ref(contact_status, data.x3d_ref);
    data.diff_3d = robot.framePosition(contact_frame_id_) - data.x3d_ref;
    l += (x3d_weight_.array()*data.diff_3d.array()*data.diff_3d.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0.0;
  }
}


void SwingFootCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (!contact_status.isContactActive(contact_index_)) {
    data.J_6d.setZero();
    robot.getFrameJacobian(contact_frame_id_, data.J_6d);
    data.J_3d.noalias() 
        = robot.frameRotation(contact_frame_id_) * data.J_6d.template topRows<3>();
    kkt_residual.lq().noalias() 
        += grid_info.dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.diff_3d;
  }
}


void SwingFootCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (!contact_status.isContactActive(contact_index_)) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.J_3d.transpose() * x3d_weight_.asDiagonal() * data.J_3d;
  }
}


double SwingFootCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  return 0.0;
}


void SwingFootCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  // Do nothing
}


void SwingFootCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  // Do nothing
}


double SwingFootCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  return 0.0;
}


void SwingFootCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  // Do nothing
}


void SwingFootCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  // Do nothing
}

} // namespace robotoc