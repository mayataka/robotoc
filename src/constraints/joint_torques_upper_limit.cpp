#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointTorquesUpperLimit::JointTorquesUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointEffortLimit().size()),
    dim_passive_(robot.dim_passive()),
    umax_(robot.jointEffortLimit()) {
}


JointTorquesUpperLimit::JointTorquesUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    umax_() {
}


JointTorquesUpperLimit::~JointTorquesUpperLimit() {
}


bool JointTorquesUpperLimit::isFeasible(
    const Robot& robot, ConstraintComponentData& data, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  for (int i=0; i<dimc_; ++i) {
    if (u.tail(dimc_).coeff(i) > umax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorquesUpperLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  assert(dtau > 0);
  data.slack = dtau * (umax_-u.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointTorquesUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> la, Eigen::Ref<Eigen::VectorXd> lf, 
    Eigen::Ref<Eigen::VectorXd> lq,  Eigen::Ref<Eigen::VectorXd> lv) const {
  // do nothing
}


void JointTorquesUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> lu) const {
  lu.tail(dimc_).noalias() += dtau * data.dual;
}


void JointTorquesUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, Eigen::Ref<Eigen::MatrixXd> Caa,
    Eigen::Ref<Eigen::MatrixXd> Cff, Eigen::Ref<Eigen::MatrixXd> Cqq,  
    Eigen::Ref<Eigen::MatrixXd> Cvv, Eigen::Ref<Eigen::VectorXd> la,
    Eigen::Ref<Eigen::VectorXd> lf, Eigen::Ref<Eigen::VectorXd> lq, 
    Eigen::Ref<Eigen::VectorXd> lv) const {
  // do nothing
}


void JointTorquesUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& u, Eigen::Ref<Eigen::MatrixXd> Cuu, 
    Eigen::Ref<Eigen::VectorXd> lu) const {
  for (int i=0; i<dimc_; ++i) {
    Cuu.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (u.tail(dimc_)-umax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  lu.tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointTorquesUpperLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& da, 
    const Eigen::Ref<const Eigen::VectorXd>& df, 
    const Eigen::Ref<const Eigen::VectorXd>& dq, 
    const Eigen::Ref<const Eigen::VectorXd>& dv, 
    const Eigen::Ref<const Eigen::VectorXd>& du) const {
  data.dslack = - dtau * du.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dslack, data.dual, data.duality, 
                       data.ddual);
}

double JointTorquesUpperLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (u.tail(dimc_)-umax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointTorquesUpperLimit::residualSquaredNrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (u.tail(dimc_)-umax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointTorquesUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp