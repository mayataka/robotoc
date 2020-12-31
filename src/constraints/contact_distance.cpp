#include "idocp/constraints/contact_distance.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

ContactDistance::ContactDistance(const Robot& robot, const double barrier,
                                 const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimv_(robot.dimv()),
    dimc_(robot.maxPointContacts()),
    contact_frames_(robot.contactFramesIndices()),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
}


ContactDistance::ContactDistance()
  : ConstraintComponentBase(),
    dimv_(0),
    dimc_(0),
    contact_frames_(),
    fraction_to_boundary_rate_(0) {
}


ContactDistance::~ContactDistance() {
}


bool ContactDistance::useKinematics() const {
  return true;
}


KinematicsLevel ContactDistance::kinematicsLevel() const {
  return KinematicsLevel::PositionLevel;
}


void ContactDistance::allocateExtraData(ConstraintComponentData& data) const {
  data.J.clear();
  for (int i=0; i<dimc_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(6, dimv_));
  }
}


bool ContactDistance::isFeasible(Robot& robot, ConstraintComponentData& data, 
                                 const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<dimc_; ++i) {
    if (!s.isContactActive(i)) {
      if (robot.framePosition(contact_frames_[i]).coeff(2) <= 0) {
        return false;
      }
    }
  }
  return true;
}


void ContactDistance::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<dimc_; ++i) {
    data.slack.coeffRef(i) = robot.framePosition(contact_frames_[i]).coeff(2);
  }
  setSlackAndDualPositive(data);
}


void ContactDistance::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  for (int i=0; i<dimc_; ++i) {
    if (!s.isContactActive(i)) {
      robot.getFrameJacobian(contact_frames_[i], data.J[i]);
      kkt_residual.lq().noalias() 
          -= dtau * data.dual.coeff(i) * data.J[i].row(2);
    }
  }
}


void ContactDistance::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  for (int i=0; i<dimc_; ++i) {
    if (!s.isContactActive(i)) {
      kkt_matrix.Qqq().noalias()
          += (dtau * data.dual.coeff(i) / data.slack.coeff(i))
              * data.J[i].row(2).transpose() * data.J[i].row(2);
      data.residual.coeffRef(i) 
          = - robot.framePosition(contact_frames_[i]).coeff(2) 
              + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      kkt_residual.lq().noalias()
          -= (dtau * (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
                  / data.slack.coeff(i))
              * data.J[i].row(2);
    }
  }
}


void ContactDistance::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  for (int i=0; i<dimc_; ++i) {
    if (!s.isContactActive(i)) {
      data.dslack.coeffRef(i) 
          = data.J[i].row(2).dot(d.dq()) - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
    }
    else {
      data.residual.coeffRef(i)  = 0;
      data.duality.coeffRef(i)   = 0;
      data.slack.coeffRef(i)     = 1.0;
      data.dslack.coeffRef(i)    = fraction_to_boundary_rate_;
      data.dual.coeffRef(i)      = 1.0;
      data.ddual.coeffRef(i)     = fraction_to_boundary_rate_;
    }
  }
}


void ContactDistance::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<dimc_; ++i) {
    if (!s.isContactActive(i)) {
      data.residual.coeffRef(i) 
          = - robot.framePosition(contact_frames_[i]).coeff(2) 
              + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
    }
  }
}


int ContactDistance::dimc() const {
  return dimc_;
}

} // namespace idocp