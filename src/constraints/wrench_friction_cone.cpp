#include "robotoc/constraints/wrench_friction_cone.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

WrenchFrictionCone::WrenchFrictionCone(const Robot& robot, const double mu, 
                                       const double X, const double Y)
  : ConstraintComponentBase(),
    dimv_(robot.dimv()),
    dimc_(17*robot.maxNumSurfaceContacts()),
    max_num_contacts_(robot.maxNumContacts()),
    contact_frame_(robot.contactFrames()),
    contact_types_(robot.contactTypes()),
    mu_(mu),
    X_(X),
    Y_(Y),
    cone_(Eigen::MatrixXd::Zero(17, 6)) {
  if (robot.maxNumContacts() == 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: robot.maxNumContacts() must be positive!");
  }
  if (mu <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'mu' must be positive!");
  }
  if (X <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'X' must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'Y' must be positive!");
  }
  setCone(mu, X, Y);
}


WrenchFrictionCone::WrenchFrictionCone()
  : ConstraintComponentBase(),
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


WrenchFrictionCone::~WrenchFrictionCone() {
}


void WrenchFrictionCone::setFrictionCoefficient(const double mu) {
  if (mu <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'mu' must be positive!");
  }
  mu_ = mu;
  setCone(mu_, X_, Y_);
}


void WrenchFrictionCone::setRectangular(const double X, const double Y) {
  if (X <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'X' must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range(
        "[WrenchFrictionCone] invalid argument: 'Y' must be positive!");
  }
  X_ = X;
  Y_ = Y;
  setCone(mu_, X_, Y_);
}


bool WrenchFrictionCone::useKinematics() const {
  return true;
}


KinematicsLevel WrenchFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void WrenchFrictionCone::allocateExtraData(ConstraintComponentData& data) const {
  data.r.clear();
  data.r.push_back(Eigen::VectorXd::Zero(17));
}


bool WrenchFrictionCone::isFeasible(Robot& robot, 
                                    const ContactStatus& contact_status,
                                    ConstraintComponentData& data, 
                                    const SplitSolution& s) const {
  data.residual.setZero();
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


void WrenchFrictionCone::setSlack(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  ConstraintComponentData& data, 
                                  const SplitSolution& s) const {
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


void WrenchFrictionCone::evalConstraint(Robot& robot, 
                                        const ContactStatus& contact_status, 
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  data.residual.setZero();
  data.cmpl.setZero();
  data.log_barrier = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


void WrenchFrictionCone::evalDerivatives(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         ConstraintComponentData& data, 
                                         const SplitSolution& s, 
                                         SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


void WrenchFrictionCone::condenseSlackAndDual(const ContactStatus& contact_status, 
                                              ConstraintComponentData& data, 
                                              SplitKKTMatrix& kkt_matrix, 
                                              SplitKKTResidual& kkt_residual) const {
  data.cond.setZero();
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<max_num_contacts_; ++i) {
    switch (contact_types_[i]) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


void WrenchFrictionCone::expandSlackAndDual(const ContactStatus& contact_status,
                                            ConstraintComponentData& data, 
                                            const SplitDirection& d) const {
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
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
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


int WrenchFrictionCone::dimc() const {
  return dimc_;
}


void WrenchFrictionCone::setCone(const double mu, const double X, 
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