#include "idocp/constraints/joint_position_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointPositionUpperLimit::JointPositionUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.lowerJointPositionLimit().size()),
    dim_passive_(robot.dim_passive()),
    qmax_(robot.upperJointPositionLimit()) {
}


JointPositionUpperLimit::JointPositionUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    qmax_() {
}


JointPositionUpperLimit::~JointPositionUpperLimit() {
}


bool JointPositionUpperLimit::isFeasible(
    const Robot& robot, ConstraintComponentData& data, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  for (int i=0; i<dimc_; ++i) {
    if (q.tail(dimc_).coeff(i) > qmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointPositionUpperLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  assert(dtau > 0);
  data.slack = dtau * (qmax_-q.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointPositionUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> la, Eigen::Ref<Eigen::VectorXd> lf, 
    Eigen::Ref<Eigen::VectorXd> lq,  Eigen::Ref<Eigen::VectorXd> lv) const {
  lq.tail(dimc_).noalias() += dtau * data.dual;
}


void JointPositionUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> lu) const {
  // do nothing
}


void JointPositionUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, Eigen::Ref<Eigen::MatrixXd> Caa,
    Eigen::Ref<Eigen::MatrixXd> Cff, Eigen::Ref<Eigen::MatrixXd> Cqq,  
    Eigen::Ref<Eigen::MatrixXd> Cvv, Eigen::Ref<Eigen::VectorXd> la,
    Eigen::Ref<Eigen::VectorXd> lf, Eigen::Ref<Eigen::VectorXd> lq, 
    Eigen::Ref<Eigen::VectorXd> lv) const {
  for (int i=0; i<dimc_; ++i) {
    Cqq.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (q.tail(dimc_)-qmax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  lq.tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointPositionUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& u, Eigen::Ref<Eigen::MatrixXd> Cuu, 
    Eigen::Ref<Eigen::VectorXd> Cu) const {
  // do nothing
}


void JointPositionUpperLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& da, 
    const Eigen::Ref<const Eigen::VectorXd>& df, 
    const Eigen::Ref<const Eigen::VectorXd>& dq, 
    const Eigen::Ref<const Eigen::VectorXd>& dv, 
    const Eigen::Ref<const Eigen::VectorXd>& du) const {
  data.dslack = - dtau * dq.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dslack, data.dual, data.duality, 
                       data.ddual);
}

double JointPositionUpperLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (q.tail(dimc_)-qmax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointPositionUpperLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (q.tail(dimc_)-qmax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointPositionUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp