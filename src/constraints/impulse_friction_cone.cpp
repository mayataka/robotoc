#include "idocp/constraints/impulse_friction_cone.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

ImpulseFrictionCone::ImpulseFrictionCone(const Robot& robot, const double mu,
                                         const double barrier,
                                         const double fraction_to_boundary_rate)
  : ImpulseConstraintComponentBase(barrier, fraction_to_boundary_rate),
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


ImpulseFrictionCone::ImpulseFrictionCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    max_point_contacts_(0),
    dimv_(0),
    contact_frame_(),
    mu_(0),
    cone_() {
}


ImpulseFrictionCone::~ImpulseFrictionCone() {
}


void ImpulseFrictionCone::setFrictionCoefficient(const double mu) {
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


KinematicsLevel ImpulseFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void ImpulseFrictionCone::allocateExtraData(
    ConstraintComponentData& data) const {
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


bool ImpulseFrictionCone::isFeasible(Robot& robot, 
                                     ConstraintComponentData& data, 
                                     const ImpulseSplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
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


void ImpulseFrictionCone::setSlack(Robot& robot, ConstraintComponentData& data, 
                                   const ImpulseSplitSolution& s) const {
  robot.updateFrameKinematics(s.q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    Eigen::VectorXd& fWi = fW(data, i);
    fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
    frictionConeResidual(mu_, fWi, data.residual.template segment<5>(5*i));
    data.slack.template segment<5>(5*i)
        = - data.residual.template segment<5>(5*i);
  }
}


void ImpulseFrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  data.residual.setZero();
  data.cmpl.setZero();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx = 5*i;
      // Contact force expressed in the world frame.
      Eigen::VectorXd& fWi = fW(data, i);
      fLocal2World(robot, contact_frame_[i], s.f[i], fWi);
      frictionConeResidual(mu_, fWi, data.residual.template segment<5>(5*i));
      data.residual.template segment<5>(idx).noalias()
          += data.slack.template segment<5>(idx);
      computeComplementarySlackness(data, idx, 5);
    }
  }
}


void ImpulseFrictionCone::computePrimalResidualDerivatives(
    Robot& robot, ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      // Contact force expressed in the world frame.
      const Eigen::VectorXd& fWi = fW(data, i);
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
          += dgi_dq.transpose() * data.dual.template segment<5>(5*i);
      // Jacobian of the frition cone constraint with respect to the contact
      // force expressed in the local frame.
      Eigen::MatrixXd& dgi_df = dg_df(data, i);
      dgi_df.noalias() = cone_ * robot.frameRotation(contact_frame_[i]);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dgi_df.transpose() * data.dual.template segment<5>(5*i);
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx = 5*i;
      const Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      const Eigen::MatrixXd& dgi_df = dg_df(data, i);
      Eigen::VectorXd& ri = r(data, i);
      ri.array() = (data.dual.template segment<5>(idx).array()
                    *data.residual.template segment<5>(idx).array()
                    -data.cmpl.template segment<5>(idx).array())
                    / data.slack.template segment<5>(idx).array();
      kkt_residual.lq().noalias() += dgi_dq.transpose() * ri;
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dgi_df.transpose() * ri;
      Eigen::MatrixXd& dfWi_dq = dfW_dq(data, i);
      Eigen::MatrixXd& r_dgi_df = r_dg_df(data, i);
      ri.array() = data.dual.template segment<5>(idx).array() 
                    / data.slack.template segment<5>(idx).array();
      dfWi_dq.template topRows<5>().noalias() = ri.asDiagonal() * dgi_dq;
      r_dgi_df.noalias() = ri.asDiagonal() * dgi_df;
      kkt_matrix.Qqq().noalias()
          += dgi_dq.transpose() * dfWi_dq.template topRows<5>();
      kkt_matrix.Qqf().template middleCols<3>(dimf_stack).noalias()
          += dgi_dq.transpose() * r_dgi_df; 
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dgi_df.transpose() * r_dgi_df;
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::expandSlackAndDual(
    ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    const ImpulseSplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // positive data.dslack and data.ddual do not affect the step size 
  // determined by the fraction-to-boundary-rule.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<max_point_contacts_; ++i) {
    if (s.isImpulseActive(i)) {
      const int idx = 5*i;
      const Eigen::MatrixXd& dgi_dq = dg_dq(data, i);
      const Eigen::MatrixXd& dgi_df = dg_df(data, i);
      data.dslack.template segment<5>(idx).noalias()
          = - dgi_dq * d.dq() - dgi_df * d.df().template segment<3>(dimf_stack) 
            - data.residual.template segment<5>(idx);
      computeDualDirection(data, idx, 5);
      dimf_stack += 3;
    }
  }
}


int ImpulseFrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp