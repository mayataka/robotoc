#include "idocp/cost/trotting_foot_height_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TrottingFootHeightCost::TrottingFootHeightCost(const Robot& robot, 
                                               const int knee_joint_id, 
                                               const int hip_joint_id, 
                                               const double foot_length, 
                                               const double thigh_length)
  : CostFunctionComponentBase(),
    knee_joint_id_(knee_joint_id), 
    hip_joint_id_(hip_joint_id),
    foot_length_(foot_length),
    thigh_length_(thigh_length),
    stance_height_(0),
    swing_height_(0),
    t_start_(0),
    t_period_(0),
    q_weight_(0),
    qf_weight_(0),
    qi_weight_(0) {
}


TrottingFootHeightCost::TrottingFootHeightCost()
  : CostFunctionComponentBase(),
    knee_joint_id_(0), 
    hip_joint_id_(0),
    foot_length_(0),
    thigh_length_(0),
    stance_height_(0),
    swing_height_(0),
    t_start_(0),
    t_period_(0),
    q_weight_(0),
    qf_weight_(0),
    qi_weight_(0) {
}


TrottingFootHeightCost::~TrottingFootHeightCost() {
}


bool TrottingFootHeightCost::useKinematics() const {
  return false;
}


void TrottingFootHeightCost::set_ref(const Eigen::VectorXd& q_standing, 
                                     const double swing_height) {
  try {
    if (swing_height <= 0) {
      throw std::runtime_error(
          "invalid argument: swing_height must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  stance_height_ = dist(q_standing);
  swing_height_ = swing_height;
}


void TrottingFootHeightCost::set_period(const double t_start, 
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


void TrottingFootHeightCost::set_q_weight(const double q_weight) {
  q_weight_ = q_weight;
}


void TrottingFootHeightCost::set_qf_weight(const double qf_weight) {
  qf_weight_ = qf_weight;
}


void TrottingFootHeightCost::set_qi_weight(const double qi_weight) {
  qi_weight_ = qi_weight;
}


double TrottingFootHeightCost::computeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s) const {
  const double diff = hdiff(t, s.q);
  return 0.5 * dtau * q_weight_ * diff * diff;
}


double TrottingFootHeightCost::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  const double diff = hdiff(t, s.q);
  return 0.5 * qf_weight_ * diff * diff;
}


double TrottingFootHeightCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  const double diff = hdiff(t, s.q);
  return 0.5 * qi_weight_ * diff * diff;
}


void TrottingFootHeightCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  const double diff = hdiff(t, s.q);
  kkt_residual.lq().coeffRef(hip_joint_id_) 
      += dtau * q_weight_ * diff * ddist_dth1(t, s.q);
  kkt_residual.lq().coeffRef(knee_joint_id_) 
      += dtau * q_weight_ * diff * ddist_dth2(t, s.q);
}


void TrottingFootHeightCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  const double diff = hdiff(t, s.q);
  kkt_residual.lq().coeffRef(hip_joint_id_) 
      += qf_weight_ * diff * ddist_dth1(t, s.q);
  kkt_residual.lq().coeffRef(knee_joint_id_) 
      += qf_weight_ * diff * ddist_dth2(t, s.q);
}


void TrottingFootHeightCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  const double diff = hdiff(t, s.q);
  kkt_residual.lq().coeffRef(hip_joint_id_) 
      += qi_weight_ * diff * ddist_dth1(t, s.q);
  kkt_residual.lq().coeffRef(knee_joint_id_) 
      += qi_weight_ * diff * ddist_dth2(t, s.q);
}


void TrottingFootHeightCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  const double diff = hdiff(t, s.q);
  const double ddist1 = ddist_dth1(t, s.q);
  const double ddist2 = ddist_dth2(t, s.q);
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, hip_joint_id_)
      += dtau * q_weight_ * (ddist1 * ddist1 + diff * dddist_dth11(t, s.q));
  const double Qqq12 
      = dtau * q_weight_ * (ddist1 * ddist2 + diff * dddist_dth12(t, s.q));
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, knee_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , hip_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , knee_joint_id_)
      += dtau * q_weight_ * (ddist2 * ddist2 + diff * dddist_dth22(t, s.q));
}


void TrottingFootHeightCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  const double diff = hdiff(t, s.q);
  const double ddist1 = ddist_dth1(t, s.q);
  const double ddist2 = ddist_dth2(t, s.q);
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, hip_joint_id_)
      += qf_weight_ * (ddist1 * ddist1 + diff * dddist_dth11(t, s.q));
  const double Qqq12 
      = qf_weight_ * (ddist1 * ddist2 + diff * dddist_dth12(t, s.q));
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, knee_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , hip_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , knee_joint_id_)
      += qf_weight_ * (ddist2 * ddist2 + diff * dddist_dth22(t, s.q));
}


void TrottingFootHeightCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  const double diff = hdiff(t, s.q);
  const double ddist1 = ddist_dth1(t, s.q);
  const double ddist2 = ddist_dth2(t, s.q);
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, hip_joint_id_)
      += qi_weight_ * (ddist1 * ddist1 + diff * dddist_dth11(t, s.q));
  const double Qqq12 
      = qi_weight_ * (ddist1 * ddist2 + diff * dddist_dth12(t, s.q));
  kkt_matrix.Qqq().coeffRef(hip_joint_id_, knee_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , hip_joint_id_)
      += Qqq12;
  kkt_matrix.Qqq().coeffRef(knee_joint_id_ , knee_joint_id_)
      += qi_weight_ * (ddist2 * ddist2 + diff * dddist_dth22(t, s.q));
}

} // namespace idocp