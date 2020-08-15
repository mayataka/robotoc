#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include <assert.h>


namespace idocp {

JointVelocityUpperLimit::JointVelocityUpperLimit(
    const Robot& robot, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.jointVelocityLimit().size()),
    dim_passive_(robot.dim_passive()),
    vmax_(robot.jointVelocityLimit()) {
}


JointVelocityUpperLimit::JointVelocityUpperLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    vmax_() {
}


JointVelocityUpperLimit::~JointVelocityUpperLimit() {
}


bool JointVelocityUpperLimit::isFeasible(
    const Robot& robot, ConstraintComponentData& data, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  for (int i=0; i<dimc_; ++i) {
    if (v.tail(dimc_).coeff(i) > vmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityUpperLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  assert(dtau > 0);
  data.slack = dtau * (vmax_-v.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointVelocityUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> la, Eigen::Ref<Eigen::VectorXd> lf, 
    Eigen::Ref<Eigen::VectorXd> lq,  Eigen::Ref<Eigen::VectorXd> lv) const {
  lv.tail(dimc_).noalias() += dtau * data.dual;
}


void JointVelocityUpperLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> lu) const {
  // do nothing
}


void JointVelocityUpperLimit::condenseSlackAndDual(
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
    Cvv.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (v.tail(dimc_)-vmax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  lv.tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointVelocityUpperLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& u, Eigen::Ref<Eigen::MatrixXd> Cuu, 
    Eigen::Ref<Eigen::VectorXd> Cu) const {
  // do nothing
}


void JointVelocityUpperLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& da, 
    const Eigen::Ref<const Eigen::VectorXd>& df, 
    const Eigen::Ref<const Eigen::VectorXd>& dq, 
    const Eigen::Ref<const Eigen::VectorXd>& dv, 
    const Eigen::Ref<const Eigen::VectorXd>& du) const {
  data.dslack = - dtau * dv.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dslack, data.dual, data.duality, 
                       data.ddual);
}

double JointVelocityUpperLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (v.tail(dimc_)-vmax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointVelocityUpperLimit::residualSquaredNrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (v.tail(dimc_)-vmax_) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointVelocityUpperLimit::dimc() const {
  return dimc_;
}

} // namespace idocp