#include "robotoc/cost/time_varying_configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace robotoc {

TimeVaryingConfigurationSpaceCost::TimeVaryingConfigurationSpaceCost(
    const Robot& robot,
    const std::shared_ptr<TimeVaryingConfigurationRefBase>& q_ref) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    q_ref_(q_ref),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qi_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


TimeVaryingConfigurationSpaceCost::TimeVaryingConfigurationSpaceCost()
  : dimq_(0),
    dimv_(0),
    q_ref_(),
    q_weight_(),
    qf_weight_(),
    qi_weight_() {
}


TimeVaryingConfigurationSpaceCost::~TimeVaryingConfigurationSpaceCost() {
}


void TimeVaryingConfigurationSpaceCost::set_q_ref(
    const std::shared_ptr<TimeVaryingConfigurationRefBase>& q_ref) {
  q_ref_ = q_ref;
}


void TimeVaryingConfigurationSpaceCost::set_q_weight(
    const Eigen::VectorXd& q_weight) {
  try {
    if (q_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  q_weight_ = q_weight;
}


void TimeVaryingConfigurationSpaceCost::set_q_weight_terminal(
    const Eigen::VectorXd& q_weight_terminal) {
  try {
    if (q_weight_terminal.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight_terminal.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qf_weight_ = q_weight_terminal;
}


void TimeVaryingConfigurationSpaceCost::set_q_weight_impulse(
    const Eigen::VectorXd& q_weight_impulse) {
  try {
    if (q_weight_impulse.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight_impulse.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qi_weight_ = q_weight_impulse;
}


bool TimeVaryingConfigurationSpaceCost::useKinematics() const {
  return false;
}


double TimeVaryingConfigurationSpaceCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s) const {
  if (q_ref_->isActive(grid_info)) {
    q_ref_->update_q_ref(robot, grid_info, data.q_ref);
    double l = 0;
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (q_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
    return 0.5 * grid_info.dt * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingConfigurationSpaceCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  if (q_ref_->isActive(grid_info)) {
    q_ref_->update_q_ref(robot, grid_info, data.q_ref);
    if (robot.hasFloatingBase()) {
      robot.dSubtractConfiguration_dqf(s.q, data.q_ref, data.J_qdiff);
      kkt_residual.lq().noalias()
          += grid_info.dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += grid_info.dt * q_weight_.array() * data.qdiff.array();
    }
  }
}


void TimeVaryingConfigurationSpaceCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  if (q_ref_->isActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias()
          += grid_info.dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += grid_info.dt * q_weight_;
    }
  }
}


double TimeVaryingConfigurationSpaceCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  if (q_ref_->isActive(grid_info)) {
    q_ref_->update_q_ref(robot, grid_info, data.q_ref);
    double l = 0;
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qf_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingConfigurationSpaceCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (q_ref_->isActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      robot.dSubtractConfiguration_dqf(s.q, data.q_ref, data.J_qdiff);
      kkt_residual.lq().noalias()
          += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += qf_weight_.array() * data.qdiff.array();
    }
  }
}


void TimeVaryingConfigurationSpaceCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (q_ref_->isActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias()
          += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
    }
  }
}


double TimeVaryingConfigurationSpaceCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  if (q_ref_->isActive(grid_info)) {
    q_ref_->update_q_ref(robot, grid_info, data.q_ref);
    double l = 0;
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qi_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
    return 0.5 * l;
  }
  else {
    return 0;
  }
}


void TimeVaryingConfigurationSpaceCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (q_ref_->isActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      robot.dSubtractConfiguration_dqf(s.q, data.q_ref, data.J_qdiff);
      kkt_residual.lq().noalias()
          += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += qi_weight_.array() * data.qdiff.array();
    }
  }
}


void TimeVaryingConfigurationSpaceCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (q_ref_->isActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias()
          += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += qi_weight_;
    }
  }
}

} // namespace robotoc