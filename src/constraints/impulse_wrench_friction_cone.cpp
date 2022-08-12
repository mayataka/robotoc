#include "robotoc/constraints/impulse_wrench_friction_cone.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

ImpulseWrenchFrictionCone::ImpulseWrenchFrictionCone(
    const Robot& robot, const double mu, const double X, const double Y)
  : ImpulseConstraintComponentBase(),
    dimc_(17*robot.maxNumSurfaceContacts()),
    max_num_contacts_(robot.maxNumContacts()),
    contact_frame_(robot.contactFrames()),
    contact_types_(robot.contactTypes()),
    mu_(mu),
    X_(X),
    Y_(Y),
    cone_(Eigen::MatrixXd::Zero(17, 6)) {
  if (mu <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: mu must be positive!");
  }
  if (X <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: X must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: Y must be positive!");
  }
  setCone(mu, X, Y);
}


ImpulseWrenchFrictionCone::ImpulseWrenchFrictionCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    max_num_contacts_(0),
    dimv_(0),
    contact_frame_(),
    contact_types_(),
    mu_(0),
    X_(0),
    Y_(0),
    cone_() {
}


ImpulseWrenchFrictionCone::~ImpulseWrenchFrictionCone() {
}


void ImpulseWrenchFrictionCone::setFrictionCoefficient(const double mu) {
  if (mu <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: mu must be positive!");
  }
  mu_ = mu;
  cone_ <<  0,  0, -1, 
            1,  0, -(mu/std::sqrt(2)),
           -1,  0, -(mu/std::sqrt(2)),
            0,  1, -(mu/std::sqrt(2)),
            0, -1, -(mu/std::sqrt(2));
}


void ImpulseWrenchFrictionCone::setRectangular(const double X, const double Y) {
  if (X <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: X must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range("[ImpulseWrenchFrictionCone] invalid argument: Y must be positive!");
  }
  X_ = X;
  Y_ = Y;
  setCone(mu_, X_, Y_);
}


KinematicsLevel ImpulseWrenchFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void ImpulseWrenchFrictionCone::allocateExtraData(
    ConstraintComponentData& data) const {
  data.r.clear();
  data.r.push_back(Eigen::VectorXd::Zero(17));
}


bool ImpulseWrenchFrictionCone::isFeasible(Robot& robot, 
                                           const ImpulseStatus& impulse_status,
                                           ConstraintComponentData& data, 
                                           const ImpulseSplitSolution& s) const {
  data.residual.setZero();
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data.residual.template segment<17>(c_begin).noalias() 
              += cone_ * s.f[i];
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  if (data.residual.maxCoeff() > 0) {
    return false;
  }
  else {
    return true;
  }
}


void ImpulseWrenchFrictionCone::setSlack(Robot& robot, 
                                         const ImpulseStatus& impulse_status,
                                         ConstraintComponentData& data, 
                                         const ImpulseSplitSolution& s) const {
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data.residual.template segment<17>(c_begin).noalias() = cone_ * s.f[i];
          data.slack.template segment<17>(c_begin)
              = - data.residual.template segment<17>(c_begin);
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchFrictionCone::evalConstraint(Robot& robot, 
                                               const ImpulseStatus& impulse_status,
                                               ConstraintComponentData& data, 
                                               const ImpulseSplitSolution& s) const {
  data.residual.setZero();
  data.cmpl.setZero();
  data.log_barrier = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data.residual.template segment<17>(c_begin).noalias() 
              = cone_ * s.f[i] + data.slack.template segment<17>(c_begin);
          computeComplementarySlackness<17>(data, c_begin);
          data.log_barrier += logBarrier(data.slack.template segment<17>(c_begin));
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchFrictionCone::evalDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status,
    ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          kkt_residual.lf().template segment<6>(dimf_stack).noalias()
              += cone_.transpose() * data.dual.template segment<17>(c_begin);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchFrictionCone::condenseSlackAndDual(
    const ImpulseStatus& impulse_status, ConstraintComponentData& data, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  data.cond.setZero();
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data.r[0].array() = data.dual.template segment<17>(c_begin).array() 
                                / data.slack.template segment<17>(c_begin).array();
          kkt_matrix.Qff().template block<6, 6>(dimf_stack, dimf_stack).noalias()
              += cone_.transpose() * data.r[0].asDiagonal() * cone_;
          computeCondensingCoeffcient<17>(data, c_begin);
          kkt_residual.lf().template segment<6>(dimf_stack).noalias()
              += cone_.transpose() * data.cond.template segment<17>(c_begin);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchFrictionCone::expandSlackAndDual(
     const ImpulseStatus& impulse_status, ConstraintComponentData& data, 
    const ImpulseSplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // positive data.dslack and data.ddual do not affect the step size 
  // determined by the fraction-to-boundary-rule.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int c_begin = 0;
  int dimf_stack = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data.dslack.template segment<17>(c_begin).noalias()
              = - cone_ * d.df().template segment<6>(dimf_stack) 
                - data.residual.template segment<17>(c_begin);
          computeDualDirection<17>(data, c_begin);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


int ImpulseWrenchFrictionCone::dimc() const {
  return dimc_;
}


void ImpulseWrenchFrictionCone::setCone(const double mu, const double X, 
                                        const double Y) {
  cone_ <<  0,  0, -1, 0, 0, 0,
           -1,  0, -mu, 0, 0, 0,
            1,  0, -mu, 0, 0, 0,
            0, -1, -mu, 0, 0, 0,
            0,  1, -mu, 0, 0, 0,
            0,  0, -Y, -1,  0, 0,
            0,  0, -Y,  1,  0, 0,
            0,  0, -X,  0, -1, 0,
            0,  0, -X,  0,  1, 0,
            -Y, -X, -(X+Y)*mu,  mu,  mu, -1,
            -Y,  X, -(X+Y)*mu,  mu, -mu, -1,
             Y, -X, -(X+Y)*mu, -mu,  mu, -1,
             Y,  X, -(X+Y)*mu, -mu, -mu, -1,
             Y,  X, -(X+Y)*mu,  mu,  mu,  1,
             Y, -X, -(X+Y)*mu,  mu, -mu,  1,
            -Y,  X, -(X+Y)*mu, -mu,  mu,  1,
            -Y, -X, -(X+Y)*mu, -mu, -mu,  1;
}

} // namespace robotoc