#include "robotoc/constraints/impulse_wrench_cone.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

ImpulseWrenchCone::ImpulseWrenchCone(const Robot& robot, 
                                     const double X, const double Y)
  : ImpulseConstraintComponentBase(),
    dimc_(17*robot.maxNumSurfaceContacts()),
    max_num_contacts_(robot.maxNumContacts()),
    contact_frame_(robot.contactFrames()),
    contact_types_(robot.contactTypes()),
    X_(X),
    Y_(Y) {
  if (robot.maxNumContacts() == 0) {
    throw std::out_of_range(
        "[ImpulseWrenchCone] invalid argument: robot.maxNumContacts() must be positive!");
  }
  if (X <= 0) {
    throw std::out_of_range(
        "[ImpulseWrenchCone] invalid argument: 'X' must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range(
        "[ImpulseWrenchCone] invalid argument: 'Y' must be positive!");
  }
}


ImpulseWrenchCone::ImpulseWrenchCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    max_num_contacts_(0),
    dimv_(0),
    contact_frame_(),
    contact_types_(),
    X_(0),
    Y_(0) {
}


ImpulseWrenchCone::~ImpulseWrenchCone() {
}


void ImpulseWrenchCone::setRectangular(const double X, const double Y) {
  if (X <= 0) {
    throw std::out_of_range(
        "[ImpulseWrenchCone] invalid argument: 'X' must be positive!");
  }
  if (Y <= 0) {
    throw std::out_of_range(
        "[ImpulseWrenchCone] invalid argument: 'Y' must be positive!");
  }
  X_ = X;
  Y_ = Y;
}


KinematicsLevel ImpulseWrenchCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void ImpulseWrenchCone::allocateExtraData(ConstraintComponentData& data) const {
  data.r.clear();
  data.r.push_back(Eigen::VectorXd::Zero(17));
  const double mu = 0.7;
  Eigen::MatrixXd cone = Eigen::MatrixXd::Zero(17, 6);
  computeCone(mu, cone);
  data.J.clear();
  for (int i=0; i<max_num_contacts_; ++i) {
    data.J.push_back(cone);
  }
}


bool ImpulseWrenchCone::isFeasible(Robot& robot, 
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
          Eigen::MatrixXd& cone_i = data.J[i];
          updateCone(impulse_status.frictionCoefficient(i), cone_i);
          data.residual.template segment<17>(c_begin).noalias() 
              += cone_i * s.f[i];
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


void ImpulseWrenchCone::setSlack(Robot& robot, 
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
          Eigen::MatrixXd& cone_i = data.J[i];
          updateCone(impulse_status.frictionCoefficient(i), cone_i);
          data.residual.template segment<17>(c_begin).noalias() = cone_i * s.f[i];
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


void ImpulseWrenchCone::evalConstraint(Robot& robot, 
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
          Eigen::MatrixXd& cone_i = data.J[i];
          updateCone(impulse_status.frictionCoefficient(i), cone_i);
          data.residual.template segment<17>(c_begin).noalias() 
              = cone_i * s.f[i] + data.slack.template segment<17>(c_begin);
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


void ImpulseWrenchCone::evalDerivatives(Robot& robot, 
                                        const ImpulseStatus& impulse_status, 
                                        ConstraintComponentData& data, 
                                        const ImpulseSplitSolution& s, 
                                        SplitKKTResidual& kkt_residual) const {
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
          const Eigen::MatrixXd& cone_i = data.J[i];
          kkt_residual.lf().template segment<6>(dimf_stack).noalias()
              += cone_i.transpose() * data.dual.template segment<17>(c_begin);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchCone::condenseSlackAndDual(const ImpulseStatus& impulse_status, 
                                             ConstraintComponentData& data, 
                                             SplitKKTMatrix& kkt_matrix, 
                                             SplitKKTResidual& kkt_residual) const {
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
          const Eigen::MatrixXd& cone_i = data.J[i];
          data.r[0].array() = data.dual.template segment<17>(c_begin).array() 
                                / data.slack.template segment<17>(c_begin).array();
          kkt_matrix.Qff().template block<6, 6>(dimf_stack, dimf_stack).noalias()
              += cone_i.transpose() * data.r[0].asDiagonal() * cone_i;
          computeCondensingCoeffcient<17>(data, c_begin);
          kkt_residual.lf().template segment<6>(dimf_stack).noalias()
              += cone_i.transpose() * data.cond.template segment<17>(c_begin);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
}


void ImpulseWrenchCone::expandSlackAndDual(const ImpulseStatus& impulse_status, 
                                           ConstraintComponentData& data, 
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
          const Eigen::MatrixXd& cone_i = data.J[i];
          data.dslack.template segment<17>(c_begin).noalias()
              = - cone_i * d.df().template segment<6>(dimf_stack) 
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


int ImpulseWrenchCone::dimc() const {
  return dimc_;
}


void ImpulseWrenchCone::computeCone(const double mu, Eigen::MatrixXd& cone) const {
  const double XYmu = (X_+Y_)*mu;
  cone.resize(17, 6);
  cone <<  0,  0,  -1,  0,  0,  0,
          -1,  0, -mu,  0,  0,  0,
           1,  0, -mu,  0,  0,  0,
           0, -1, -mu,  0,  0,  0,
           0,  1, -mu,  0,  0,  0,
           0,  0, -Y_, -1,  0,  0,
           0,  0, -Y_,  1,  0,  0,
           0,  0, -X_,  0, -1,  0,
           0,  0, -X_,  0,  1,  0,
          -Y_, -X_, -XYmu,  mu,  mu, -1,
          -Y_,  X_, -XYmu,  mu, -mu, -1,
           Y_, -X_, -XYmu, -mu,  mu, -1,
           Y_,  X_, -XYmu, -mu, -mu, -1,
           Y_,  X_, -XYmu,  mu,  mu,  1,
           Y_, -X_, -XYmu,  mu, -mu,  1,
          -Y_,  X_, -XYmu, -mu,  mu,  1,
          -Y_, -X_, -XYmu, -mu, -mu,  1;
}


void ImpulseWrenchCone::updateCone(const double mu, Eigen::MatrixXd& cone) const {
  for (int i=1; i<5; ++i) {
    cone.coeffRef(i, 2) = -mu;
  }
  cone.coeffRef(5, 2) = -Y_;
  cone.coeffRef(6, 2) = -Y_;
  cone.coeffRef(7, 2) = -X_;
  cone.coeffRef(8, 2) = -X_;
  const double XYmu = (X_+Y_)*mu;
  cone.row(9)  << -Y_, -X_, -XYmu,  mu,  mu, -1;
  cone.row(10) << -Y_,  X_, -XYmu,  mu, -mu, -1;
  cone.row(11) <<  Y_, -X_, -XYmu, -mu,  mu, -1;
  cone.row(12) <<  Y_,  X_, -XYmu, -mu, -mu, -1;
  cone.row(13) <<  Y_,  X_, -XYmu,  mu,  mu,  1;
  cone.row(14) <<  Y_, -X_, -XYmu,  mu, -mu,  1;
  cone.row(15) << -Y_,  X_, -XYmu, -mu,  mu,  1;
  cone.row(16) << -Y_, -X_, -XYmu, -mu, -mu,  1;
}

} // namespace robotoc