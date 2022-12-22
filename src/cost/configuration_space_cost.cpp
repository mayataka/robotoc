#include "robotoc/cost/configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace robotoc {

ConfigurationSpaceCost::ConfigurationSpaceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    q_ref_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    u_ref_(Eigen::VectorXd::Zero(robot.dimu())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    u_weight_(Eigen::VectorXd::Zero(robot.dimu())),
    q_weight_terminal_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_terminal_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_impact_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_impact_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_weight_impact_(Eigen::VectorXd::Zero(robot.dimv())),
    ref_(),
    use_nonconst_ref_(false),
    enable_q_cost_(false), 
    enable_v_cost_(false), 
    enable_a_cost_(false), 
    enable_u_cost_(false), 
    enable_q_cost_terminal_(false), 
    enable_v_cost_terminal_(false),
    enable_q_cost_impact_(false), 
    enable_v_cost_impact_(false), 
    enable_dv_cost_impact_(false) {
  if (robot.hasFloatingBase()) {
    robot.normalizeConfiguration(q_ref_);
  }
}


ConfigurationSpaceCost::ConfigurationSpaceCost(
    const Robot& robot, const std::shared_ptr<ConfigurationSpaceRefBase>& ref) 
  : ConfigurationSpaceCost(robot) {
  set_ref(ref);
}


ConfigurationSpaceCost::ConfigurationSpaceCost()
  : CostFunctionComponentBase(),
    dimq_(0),
    dimv_(0),
    dimu_(0),
    q_ref_(),
    v_ref_(),
    u_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    u_weight_(),
    q_weight_terminal_(),
    v_weight_terminal_(),
    q_weight_impact_(),
    v_weight_impact_(),
    dv_weight_impact_(),
    ref_(),
    use_nonconst_ref_(false),
    enable_q_cost_(false), 
    enable_v_cost_(false), 
    enable_a_cost_(false), 
    enable_u_cost_(false), 
    enable_q_cost_terminal_(false), 
    enable_v_cost_terminal_(false),
    enable_q_cost_impact_(false), 
    enable_v_cost_impact_(false), 
    enable_dv_cost_impact_(false) {
}


ConfigurationSpaceCost::~ConfigurationSpaceCost() {
}


void ConfigurationSpaceCost::set_ref(
    const std::shared_ptr<ConfigurationSpaceRefBase>& ref) {
  ref_ = ref;
  use_nonconst_ref_ = true;
}


void ConfigurationSpaceCost::set_q_ref(const Eigen::VectorXd& q_ref) {
  if (q_ref.size() != dimq_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: q_ref.size() must be " + std::to_string(dimq_) + "!");
  }
  q_ref_ = q_ref;
  use_nonconst_ref_ = false;
}


void ConfigurationSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref) {
  if (v_ref.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: v_ref.size() must be " + std::to_string(dimv_) + "!");
  }
  v_ref_ = v_ref;
}


void ConfigurationSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref) {
  if (u_ref.size() != dimu_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: u_ref.size() must be " + std::to_string(dimu_) + "!");
  }
  u_ref_ = u_ref;
}


void ConfigurationSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight) {
  if (q_weight.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: q_weight.size() must be " + std::to_string(dimv_) + "!");
  }
  if (q_weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'q_weight' must be non-negative!");
  }
  q_weight_ = q_weight;
  enable_q_cost_ = (!q_weight.isZero());
}


void ConfigurationSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight) {
  if (v_weight.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: v_weight.size() must be " + std::to_string(dimv_) + "!");
  }
  if (v_weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'v_weight' must be non-negative!");
  }
  v_weight_ = v_weight;
  enable_v_cost_ = (!v_weight.isZero());
}


void ConfigurationSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight) {
  if (a_weight.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: a_weight.size() must be " + std::to_string(dimv_) + "!");
  }
  if (a_weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'a_weight' must be non-negative!");
  }
  a_weight_ = a_weight;
  enable_a_cost_ = (!a_weight.isZero());
}


void ConfigurationSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight) {
  if (u_weight.size() != dimu_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: u_weight.size() must be " + std::to_string(dimu_) + "!");
  }
  if (u_weight.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'u_weight' must be non-negative!");
  }
  u_weight_ = u_weight;
  enable_u_cost_ = (!u_weight.isZero());
}


void ConfigurationSpaceCost::set_q_weight_terminal(
    const Eigen::VectorXd& q_weight_terminal) {
  if (q_weight_terminal.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: q_weight_terminal.size() must be " + std::to_string(dimv_) + "!");
  }
  if (q_weight_terminal.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'q_weight_terminal' must be non-negative!");
  }
  q_weight_terminal_ = q_weight_terminal;
  enable_q_cost_terminal_ = (!q_weight_terminal.isZero());
}


void ConfigurationSpaceCost::set_v_weight_terminal(
    const Eigen::VectorXd& v_weight_terminal) {
  if (v_weight_terminal.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: v_weight_terminal.size() must be " + std::to_string(dimv_) + "!");
  }
  if (v_weight_terminal.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'v_weight_terminal' must be non-negative!");
  }
  v_weight_terminal_ = v_weight_terminal;
  enable_v_cost_terminal_ = (!v_weight_terminal.isZero());
}


void ConfigurationSpaceCost::set_q_weight_impact(
    const Eigen::VectorXd& q_weight_impact) {
  if (q_weight_impact.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: q_weight_impact.size() must be " + std::to_string(dimv_) + "!");
  }
  if (q_weight_impact.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'q_weight_impact' must be non-negative!");
  }
  q_weight_impact_ = q_weight_impact;
  enable_q_cost_impact_ = (!q_weight_impact.isZero());
}


void ConfigurationSpaceCost::set_v_weight_impact(
    const Eigen::VectorXd& v_weight_impact) {
  if (v_weight_impact.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: v_weight_impact.size() must be " + std::to_string(dimv_) + "!");
  }
  if (v_weight_impact.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'v_weight_impact' must be non-negative!");
  }
  v_weight_impact_ = v_weight_impact;
  enable_v_cost_impact_ = (!v_weight_impact.isZero());
}


void ConfigurationSpaceCost::set_dv_weight_impact(
    const Eigen::VectorXd& dv_weight_impact) {
  if (dv_weight_impact.size() != dimv_) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: dv_weight_impact.size() must be " + std::to_string(dimv_) + "!");
  }
  if (dv_weight_impact.minCoeff() < 0.0) {
    throw std::invalid_argument(
        "[ConfigurationSpaceCost] invalid argument: elements of 'dv_weight_impact' must be non-negative!");
  }
  dv_weight_impact_ = dv_weight_impact;
  enable_dv_cost_impact_ = (!dv_weight_impact.isZero());
}


double ConfigurationSpaceCost::evalStageCost(Robot& robot, 
                                             const ContactStatus& contact_status, 
                                             const GridInfo& grid_info,
                                             const SplitSolution& s,
                                             CostFunctionData& data) const {
  double l = 0;
  if (enable_q_cost_ && isCostConfigActive(grid_info)) {
    evalConfigDiff(robot, data, grid_info, s.q);
    l += (q_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  if (enable_v_cost_) {
    l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  }
  if (enable_a_cost_) {
    l += (a_weight_.array()*s.a.array()*s.a.array()).sum();
  }
  if (enable_u_cost_) {
    l += (u_weight_.array()*(s.u-u_ref_).array()*(s.u-u_ref_).array()).sum();
  }
  return 0.5 * l;
}


void ConfigurationSpaceCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTResidual& kkt_residual) const {
  if (enable_q_cost_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      evalConfigDiffJac(robot, data, grid_info, s.q);
      kkt_residual.lq().noalias()
          += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += q_weight_.array() * data.qdiff.array();
    }
  }
  if (enable_v_cost_) {
    kkt_residual.lv().array() += v_weight_.array() * (s.v.array()-v_ref_.array());
  }
  if (enable_a_cost_) {
    kkt_residual.la.array() += a_weight_.array() * s.a.array();
  }
  if (enable_u_cost_) {
    kkt_residual.lu.array() += u_weight_.array() * (s.u.array()-u_ref_.array());
  }
}


void ConfigurationSpaceCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTMatrix& kkt_matrix) const {
  if (enable_q_cost_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias() += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += q_weight_;
    }
  }
  if (enable_v_cost_) {
    kkt_matrix.Qvv().diagonal().noalias() += v_weight_;
  }
  if (enable_a_cost_) {
    kkt_matrix.Qaa.diagonal().noalias() += a_weight_;
  }
  if (enable_u_cost_) {
    kkt_matrix.Quu.diagonal().noalias() += u_weight_;
  }
}


double ConfigurationSpaceCost::evalTerminalCost(Robot& robot, 
                                                const GridInfo& grid_info, 
                                                const SplitSolution& s,
                                                CostFunctionData& data) const {
  double l = 0;
  if (enable_q_cost_terminal_ && isCostConfigActive(grid_info)) {
    evalConfigDiff(robot, data, grid_info, s.q);
    l += (q_weight_terminal_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  if (enable_v_cost_terminal_) {
    l += (v_weight_terminal_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  }
  return 0.5 * l;
}


void ConfigurationSpaceCost::evalTerminalCostDerivatives(
    Robot& robot, const GridInfo& grid_info, const SplitSolution& s, 
    CostFunctionData& data, SplitKKTResidual& kkt_residual) const {
  if (enable_q_cost_terminal_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      evalConfigDiffJac(robot, data, grid_info, s.q);
      kkt_residual.lq().noalias()
          += data.J_qdiff.transpose() * q_weight_terminal_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += q_weight_terminal_.array() * data.qdiff.array();
    }
  }
  if (enable_v_cost_terminal_) {
    kkt_residual.lv().array()
        += v_weight_terminal_.array() * (s.v.array()-v_ref_.array());
  }
}


void ConfigurationSpaceCost::evalTerminalCostHessian(
    Robot& robot, const GridInfo& grid_info, const SplitSolution& s, 
    CostFunctionData& data, SplitKKTMatrix& kkt_matrix) const {
  if (enable_q_cost_terminal_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias()
          += data.J_qdiff.transpose() * q_weight_terminal_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += q_weight_terminal_;
    }
  }
  if (enable_v_cost_terminal_) {
    kkt_matrix.Qvv().diagonal().noalias() += v_weight_terminal_;
  }
}


double ConfigurationSpaceCost::evalImpactCost(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info, 
    const SplitSolution& s, CostFunctionData& data) const {
  double l = 0;
  if (enable_q_cost_impact_ && isCostConfigActive(grid_info)) {
    evalConfigDiff(robot, data, grid_info, s.q);
    l += (q_weight_impact_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  if (enable_v_cost_impact_) {
    l += (v_weight_impact_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  }
  if (enable_dv_cost_impact_) {
    l += (dv_weight_impact_.array()*s.dv.array()*s.dv.array()).sum();
  }
  return 0.5 * l;
}


void ConfigurationSpaceCost::evalImpactCostDerivatives(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTResidual& kkt_residual) const {
  if (enable_q_cost_impact_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      evalConfigDiffJac(robot, data, grid_info, s.q);
      kkt_residual.lq().noalias()
          += data.J_qdiff.transpose() * q_weight_impact_.asDiagonal() * data.qdiff;
    }
    else {
      kkt_residual.lq().array() += q_weight_impact_.array() * data.qdiff.array();
    }
  }
  if (enable_v_cost_impact_) {
    kkt_residual.lv().array()
        += v_weight_impact_.array() * (s.v.array()-v_ref_.array());
  }
  if (enable_dv_cost_impact_) {
    kkt_residual.ldv.array() += dv_weight_impact_.array() * s.dv.array();
  }
}


void ConfigurationSpaceCost::evalImpactCostHessian(
    Robot& robot, const ImpactStatus& impact_status, const GridInfo& grid_info,
    const SplitSolution& s, CostFunctionData& data, 
    SplitKKTMatrix& kkt_matrix) const {
  if (enable_q_cost_impact_ && isCostConfigActive(grid_info)) {
    if (robot.hasFloatingBase()) {
      kkt_matrix.Qqq().noalias()
          += data.J_qdiff.transpose() * q_weight_impact_.asDiagonal() * data.J_qdiff;
    }
    else {
      kkt_matrix.Qqq().diagonal().noalias() += q_weight_impact_;
    }
  }
  if (enable_v_cost_impact_) {
    kkt_matrix.Qvv().diagonal().noalias() += v_weight_impact_;
  }
  if (enable_dv_cost_impact_) {
    kkt_matrix.Qdvdv.diagonal().noalias() += dv_weight_impact_;
  }
}

} // namespace robotoc