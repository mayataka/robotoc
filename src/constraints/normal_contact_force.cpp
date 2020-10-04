#include "idocp/constraints/normal_contact_force.hpp"

#include <exception>
#include <iostream>
#include <assert.h>


namespace idocp {

NormalCotactForce::NormalCotactForce(const Robot& robot, const double barrier, 
                                     const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.num_point_contacts()) {
}


NormalCotactForce::NormalCotactForce()
  : ConstraintComponentBase(),
    dimc_(0) {
}


NormalCotactForce::~NormalCotactForce() {
}


bool NormalCotactForce::useKinematics() const {
  return false;
}


bool NormalCotactForce::isFeasible(Robot& robot, ConstraintComponentData& data, 
                                   const SplitSolution& s) const {
  for (int i=0; i<robot.num_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      if (s.f[i].coeff(2) < 0) {
        return false;
      }
    }
  }
  return true;
}


void NormalCotactForce::setSlackAndDual(Robot& robot, 
                                        ConstraintComponentData& data, 
                                        const double dtau, 
                                        const SplitSolution& s) const {
  assert(dtau > 0);
  for (int i=0; i<robot.num_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      data.slack.coeffRef(i) = f[i].coeff(2);
    }
  }
  setSlackAndDualPositive(data);
}


void NormalCotactForce::augmentDualResidual(Robot& robot, 
                                            ConstraintComponentData& data, 
                                            const double dtau, 
                                            KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      kkt_residual.lf().coeffRef(dimf_stack+2) = - dtau * data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void NormalCotactForce::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int num_contact_stack = 0;
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      data.residual.coeffRef(num_contact_stack) 
          = data.slack.coeff(num_contact_stack) + f_tangent - mu_f_normal;
      data.duality.coeffRef(num_contact_stack)
          = data.slack.coeff(num_contact_stack) 
              * data.dual.coeff(num_contact_stack) - barrier_;

      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2)
          += dtau * dtau * data.dual.coeff(num_contact_stack) 
                            / data.slack.coeff(num_contact_stack);
      data.residual.coeffRef(num_contact_stack) 
          = - dtau * s.f[i].coeff(2) + data.slack.coeff(num_contact_stack);
      data.duality.coeffRef()

      computeDuality(data);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += df * () / 

          = 2 * dtau * s.f[i].coeff(0) * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+1) 
          = 2 * dtau * s.f[i].coeff(1) * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+2) 
          = - 2 * dtau * mu_ * mu_ * s.f[i].coeff(2) * data.dual.coeff(i);
      ++num_contact_stack;
      dimf_stack += 3;
    }
  }


  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lv().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void NormalCotactForce::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data);
}


double NormalCotactForce::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  double norm = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      norm += std::abs(slack.coeff(i) - dtau * s.f[i].coeff(2));
      norm += std::abs(slack.coeff(i) * dual.coeff(i) - barrier_);
    }
  }
  return norm;
}


double NormalCotactForce::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  double norm = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (robot.is_contact_active(i)) {
      const double residual = slack.coeff(i) - dtau * s.f[i].coeff(2);
      const double duality = slack.coeff(i) * dual.coeff(i) - barrier_;
      norm += residual * residual + duality * duality;
    }
  }
  return norm;
}


int NormalCotactForce::dimc() const {
  return dimc_;
}

} // namespace idocp