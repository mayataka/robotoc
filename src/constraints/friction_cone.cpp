#include "idocp/constraints/friction_cone.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, const double mu, 
                           const double barrier,
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimv_(robot.dimv()),
    dimc_(5*robot.maxPointContacts()),
    max_point_contacts_(robot.maxPointContacts()),
    contact_frame_(robot.contactFrames()),
    mu_(mu),
    cone_(Eigen::MatrixXd::Zero(5, 3)) {
  try {
    if (robot.maxPointContacts() == 0) {
      throw std::out_of_range(
          "Invalid argument: robot.maxPointContacts() must be positive!");
    }
    if (mu <= 0) {
      throw std::out_of_range(
          "Invalid argument: mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  cone_ <<  0,  0, -1, 
            1,  0, -(mu/std::sqrt(2)),
           -1,  0, -(mu/std::sqrt(2)),
            0,  1, -(mu/std::sqrt(2)),
            0, -1, -(mu/std::sqrt(2));
}


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    max_point_contacts_(0),
    dimv_(0),
    contact_frame_(),
    mu_(0),
    cone_() {
}


FrictionCone::~FrictionCone() {
}


void FrictionCone::setFrictionCoefficient(const double mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range("Invalid argument: mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  mu_ = mu;
  cone_ <<  0,  0, -1, 
            1,  0, -(mu/std::sqrt(2)),
           -1,  0, -(mu/std::sqrt(2)),
            0,  1, -(mu/std::sqrt(2)),
            0, -1, -(mu/std::sqrt(2));
}


bool FrictionCone::useKinematics() const {
  return true;
}


KinematicsLevel FrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void FrictionCone::allocateExtraData(ConstraintComponentData& data) const {
  data.r.clear();
  for (int i=0; i<max_point_contacts_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(3)); // fWi
  }
  for (int i=0; i<max_point_contacts_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(5)); // ri
  }
  data.J.clear();
  for (int i=0; i<max_point_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, dimv_)); // dgi_dq
  }
  for (int i=0; i<max_point_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3)); // dgi_df
  }
  for (int i=0; i<max_point_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(6, dimv_)); // dfWi_dq
  }
  for (int i=0; i<max_point_contacts_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3)); // r_dgi_df
  }
}


bool FrictionCone::isFeasible(Robot& robot, ConstraintComponentData& data, 
                              const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      Eigen::VectorXd& fWi = fW(data, i);
      fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
      frictionConeResidual(mu_, fWi, data.residual.template segment<5>(5*i));
      if (data.residual.maxCoeff() > 0) {
        return false;
      }
    }
  }
  return true;
}


void FrictionCone::setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                                   const SplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    Eigen::VectorXd& fWi = fW(data, i);
    fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
    frictionConeResidual(mu_, fWi, data.residual.template segment<5>(5*i));
    data.slack.template segment<5>(5*i)
        = - data.residual.template segment<5>(5*i);
  }
  setSlackAndDualPositive(data);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dt > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      // Contact force expressed in the world frame.
      Eigen::VectorXd& fWi = fW(data, i);
      fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
      // Jacobian of the contact force expressed in the world frame fWi 
      // with respect to the configuration q.
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      dfWi_dq.setZero();
      robot.getFrameJacobian(contact_frame_[i], dfWi_dq);
      for (int j=0; j<dimv_; ++j) {
        dfWi_dq.template topRows<3>().col(j)
            = dfWi_dq.template bottomRows<3>().col(j).cross(fWi.template head<3>());
      }
      // Jacobian of the frition cone constraint with respect to the 
      // configuration q.
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      dgi_dq.noalias() = cone_ * dfWi_dq.template topRows<3>();
      kkt_residual.lq().noalias()
          += dt * dgi_dq.transpose() * data.dual.template segment<5>(5*i);
      // Jacobian of the frition cone constraint with respect to the contact
      // force expressed in the local frame.
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
      dgi_df.noalias() = cone_ * robot.frameRotation(contact_frame_[i]);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dt * dgi_df.transpose() * data.dual.template segment<5>(5*i);
      dimf_stack += 3;
    }
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dt > 0);
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx = 5*i;
      Eigen::VectorXd& ri = r(data, i);
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
      ri.array() = (data.dual.template segment<5>(idx).array()
                    *data.residual.template segment<5>(idx).array()
                    -data.duality.template segment<5>(idx).array())
                    / data.slack.template segment<5>(idx).array();
      kkt_residual.lq().noalias() += dt * dgi_dq.transpose() * ri;
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dt * dgi_df.transpose() * ri;
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      Eigen::MatrixXd& r_dgi_df = r_dg_df(data, i);
      ri.array() = data.dual.template segment<5>(idx).array() 
                    / data.slack.template segment<5>(idx).array();
      dfWi_dq.template topRows<5>().noalias() = ri.asDiagonal() * dgi_dq;
      r_dgi_df.noalias() = ri.asDiagonal() * dgi_df;
      kkt_matrix.Qqq().noalias()
          += dt * dgi_dq.transpose() * dfWi_dq.template topRows<5>();
      kkt_matrix.Qqf().template middleCols<3>(dimf_stack).noalias()
          += dt * dgi_dq.transpose() * r_dgi_df; 
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dt * dgi_df.transpose() * r_dgi_df;
      dimf_stack += 3;
    }
  }
}


void FrictionCone::expandSlackAndDual(ConstraintComponentData& data, 
                                      const SplitSolution& s, 
                                      const SplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // - fraction_rate * (slack.coeff(i)/dslack.coeff(i)) and 
  // - fraction_rate * (dual.coeff(i)/ddual.coeff(i))  
  // at the inactive constraint index i are always negative, 
  // and therefore do not affect to step size.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isContactActive(i)) {
      const int idx = 5*i;
      Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
      data.dslack.template segment<5>(idx).noalias()
          = - dgi_dq * d.dq() - dgi_df * d.df().template segment<3>(dimf_stack) 
            - data.residual.template segment<5>(idx);
      for (int j=0; j<5; ++j) {
        data.ddual.coeffRef(idx+j) = computeDualDirection(data.slack.coeff(idx+j), 
                                                          data.dual.coeff(idx+j), 
                                                          data.dslack.coeff(idx+j), 
                                                          data.duality.coeff(idx+j));
      }
      dimf_stack += 3;
    }
  }
}


void FrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual.setZero();
  data.duality.setZero();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      Eigen::VectorXd& fWi = fW(data, i);
      fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
      frictionConeResidual(mu_, fWi, data.residual.template segment<5>(5*i));
      data.residual.template segment<5>(idx).noalias()
          += data.slack.template segment<5>(idx);
      for (int j=0; j<5; ++j) {
        data.duality.coeffRef(idx+j) = computeDuality(data.slack.coeff(idx+j), 
                                                      data.dual.coeff(idx+j));
      }
    }
  }
}


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp