#include "idocp/constraints/friction_cone.hpp"

#include <exception>
#include <iostream>
#include <assert.h>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, double mu, const double barrier, 
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.max_point_contacts()),
    mu_(mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range("invalid argment: mu must be positive");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    mu_(0) {
}


FrictionCone::~FrictionCone() {
}


bool FrictionCone::useKinematics() const {
  return false;
}


bool FrictionCone::isFeasible(Robot& robot, ConstraintComponentData& data, 
                              const SplitSolution& s) const {
  if (robot.is_contact_active(contact_index_)) {
    const double f_tangent = s.f[contact_index_].coeff(0) 
                                * s.f[contact_index_].coeff(0) 
                              + s.f[contact_index_].coeff(1) 
                                * s.f[contact_index_].coeff(1);
    const double mu_f_normal = mu_ * mu_ * s.f[contact_index_].coeff(2) 
                                         * s.f[contact_index_].coeff(2);
    if (mu_f_normal <= f_tangent) {
      return false;
    }
  }
  return true;
}


void FrictionCone::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  if (robot.is_contact_active(contact_index_)) {
    const double f_tangent = s.f[contact_index_].coeff(0) 
                                * s.f[contact_index_].coeff(0) 
                              + s.f[contact_index_].coeff(1) 
                                * s.f[contact_index_].coeff(1);
    const double mu_f_normal = mu_ * mu_ * s.f[contact_index_].coeff(2) 
                                         * s.f[contact_index_].coeff(2);
    data.slack.coeffRef(0) = mu_f_normal - f_tangent;
  }
  setSlackAndDualPositive(data.slack, data.dual);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  if (robot.is_contact_active(contact_index_)) {
    kkt_residual.lf().coeffRef(dimf_stack  ) 
        = 2 * dtau * s.f[contact_index_].coeff(0) * data.dual.coeff(0);
    kkt_residual.lf().coeffRef(dimf_stack+1) 
        = 2 * dtau * s.f[contact_index_].coeff(1) * data.dual.coeff(0);
    kkt_residual.lf().coeffRef(dimf_stack+2) 
        = - 2 * dtau * mu_ * mu_ * s.f[contact_index_].coeff(2) 
                                 * data.dual.coeff(0);
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int num_contact_stack = 0;
  int dimf_stack = 0;
  Eigen::Vector3d df;
  if (robot.is_contact_active(contact_index_)) {
    const double f_tangent = s.f[contact_index_].coeff(0) 
                                * s.f[contact_index_].coeff(0) 
                              + s.f[contact_index_].coeff(1) 
                                * s.f[contact_index_].coeff(1);
    const double mu_f_normal = mu_ * mu_ * s.f[contact_index_].coeff(2) 
                                         * s.f[contact_index_].coeff(2);

    data.residual.coeffRef(num_contact_stack) 
        = data.slack.coeff(num_contact_stack) + f_tangent - mu_f_normal;
    data.duality.coeffRef(num_contact_stack)
        = data.slack.coeff(num_contact_stack) 
            * data.dual.coeff(num_contact_stack) - barrier_;
    const double dual_per_slack =  
                                    ;
    df.coeffRef(0) = 2 * dtau * s.f[i].coeff(0);
    df.coeffRef(1) = 2 * dtau * s.f[i].coeff(1);
    df.coeffRef(2) = - 2 * dtau * mu_ * mu_ * s.f[i].coeff(2);
    kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
      += (data.dual.coeff(num_contact_stack) 
            / data.slack.coeff(num_contact_stack)) 
          * df.transpose() * df;
    kkt_residual.lf().template segment<3>(dimf_stack).noalias()
      += df * () / 

        = 2 * dtau * s.f[i].coeff(0) * data.dual.coeff(i);
    kkt_residual.lf().coeffRef(dimf_stack+1) 
        = 2 * dtau * s.f[i].coeff(1) * data.dual.coeff(i);
    kkt_residual.lf().coeffRef(dimf_stack+2) 
        = - 2 * dtau * mu_ * mu_ * s.f[i].coeff(2) * data.dual.coeff(i);
    dimf_stack += 3;
    ++num_contact_stack;
  }


  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lv().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void FrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double FrictionCone::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double FrictionCone::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp