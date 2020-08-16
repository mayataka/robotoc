#include "idocp/constraints/joint_acceleration_lower_limit.hpp"

#include <assert.h>


namespace idocp {

JointAccelerationLowerLimit::JointAccelerationLowerLimit(
    const Robot& robot, const Eigen::VectorXd& amin, const double barrier, 
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(amin.size()),
    dim_passive_(robot.dim_passive()),
    amin_(amin) {
}


JointAccelerationLowerLimit::JointAccelerationLowerLimit()
  : ConstraintComponentBase(),
    dimc_(0),
    dim_passive_(0),
    amin_() {
}


JointAccelerationLowerLimit::~JointAccelerationLowerLimit() {
}


bool JointAccelerationLowerLimit::isFeasible(
    const Robot& robot, ConstraintComponentData& data, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  for (int i=0; i<dimc_; ++i) {
    if (a.tail(dimc_).coeff(i) < amin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointAccelerationLowerLimit::setSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  assert(dtau > 0);
  data.slack = dtau * (a.tail(dimc_)-amin_);
  setSlackAndDualPositive(data.slack, data.dual);
}


void JointAccelerationLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> la, Eigen::Ref<Eigen::VectorXd> lf, 
    Eigen::Ref<Eigen::VectorXd> lq,  Eigen::Ref<Eigen::VectorXd> lv) const {
  la.tail(dimc_).noalias() -= dtau * data.dual;
}


void JointAccelerationLowerLimit::augmentDualResidual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    Eigen::Ref<Eigen::VectorXd> lu) const {
  // do nothing
}


void JointAccelerationLowerLimit::condenseSlackAndDual(
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
    Caa.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
  }
  data.residual = dtau * (amin_-a.tail(dimc_)) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  la.tail(dimc_).array() 
      -= dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void JointAccelerationLowerLimit::condenseSlackAndDual(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& u, Eigen::Ref<Eigen::MatrixXd> Cuu, 
    Eigen::Ref<Eigen::VectorXd> Cu) const {
  // do nothing
}


void JointAccelerationLowerLimit::computeSlackAndDualDirection(
    const Robot& robot, ConstraintComponentData& data, const double dtau, 
    const Eigen::Ref<const Eigen::VectorXd>& da, 
    const Eigen::Ref<const Eigen::VectorXd>& df, 
    const Eigen::Ref<const Eigen::VectorXd>& dq, 
    const Eigen::Ref<const Eigen::VectorXd>& dv, 
    const Eigen::Ref<const Eigen::VectorXd>& du) const {
  data.dslack = dtau * da.tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dslack, data.dual, data.duality, 
                       data.ddual);
}

double JointAccelerationLowerLimit::residualL1Nrom(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (amin_-a.tail(dimc_)) + data.slack;
  return data.residual.lpNorm<1>();
}


double JointAccelerationLowerLimit::squaredKKTErrorNorm(
    const Robot& robot, ConstraintComponentData& data, 
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& a, 
    const Eigen::Ref<const Eigen::VectorXd>& f, 
    const Eigen::Ref<const Eigen::VectorXd>& q, 
    const Eigen::Ref<const Eigen::VectorXd>& v, 
    const Eigen::Ref<const Eigen::VectorXd>& u) const {
  data.residual = dtau * (amin_-a.tail(dimc_)) + data.slack;
  computeDualityResidual(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int JointAccelerationLowerLimit::dimc() const {
  return dimc_;
}

} // namespace idocp