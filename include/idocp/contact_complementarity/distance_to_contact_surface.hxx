#ifndef IDOCP_DISTANCE_TO_CONTACT_SURFACE_HXX_
#define IDOCP_DISTANCE_TO_CONTACT_SURFACE_HXX_

#include "idocp/contact_complementarity/distance_to_contact_surface.hpp"

#include <exception>
#include <iostream>
#include <assert.h>


namespace idocp {

inline DistanceToContactSurface::DistanceToContactSurface(
   const Robot& robot, const double barrier, 
   const double fraction_to_boundary_rate)
  : ContactComplementarityComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.max_point_contacts()),
    frame_position_(robot.max_point_contacts(), Eigen::Vector3d::Zero()),
    frame_derivative_(robot.max_point_contacts(), 
                      Eigen::MatrixXd::Zero(3, robot.dimv())) {
}


inline DistanceToContactSurface::DistanceToContactSurface()
  : ContactComplementarityComponentBase(),
    dimc_(0),
    frame_position_(),
    frame_derivative_() {
}


inline DistanceToContactSurface::~DistanceToContactSurface() {
}


inline bool DistanceToContactSurface::isFeasible_impl(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      robot.computeContactResidual(i, frame_position_[i]);
      const double distance = frame_position_[i].coeff(2);
      if (distance < 0) {
        return false;
      }
    }
  }
  return true;
}


inline void DistanceToContactSurface::setSlackAndDual_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    robot.computeContactResidual(i, frame_position_[i]);
    const double distance = frame_position_[i].coeff(2);
    data.slack.coeffRef(i) = dtau * distance;
  }
  setSlackAndDualPositive(data.slack, data.dual);
}


inline void DistanceToContactSurface::augmentDualResidual_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      robot.computeContactDerivative(i, frame_derivative_[i]);
      kkt_residual.lq().noalias() 
          -= dtau * frame_derivative_[i].row(2) * data.dual.coeff(i);
    }
  }
}


inline void DistanceToContactSurface::condenseSlackAndDual_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      kkt_matrix.Qqq().noalias() 
          += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i) 
                         * frame_derivative_[i].row(2).transpose() 
                         * frame_derivative_[i].row(2);
      robot.computeContactResidual(i, frame_position_[i]);
      data.residual.coeffRef(i) = - dtau * frame_position_[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      kkt_residual.lq().noalias()
          -= dtau * frame_derivative_[i].row(2)
              * (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
                  / data.slack.coeff(i);
    }
  }
}


inline void DistanceToContactSurface::computeSlackAndDualDirection_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, const SplitDirection& d) const {
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      data.dslack.coeffRef(i) 
          = dtau * frame_derivative_[i].row(2).dot(d.dq()) - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
    }
  }
}


inline double DistanceToContactSurface::residualL1Nrom_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  double norm = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      robot.computeContactResidual(i, frame_position_[i]);
      const double distance = frame_position_[i].coeff(2);
      norm += std::abs(data.slack.coeff(i) - dtau * distance);
      norm += std::abs(computeDuality(data.slack.coeff(i), data.dual.coeff(i)));
    }
  }
  return norm;
}


inline double DistanceToContactSurface::squaredKKTErrorNorm_impl(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  double norm = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      robot.computeContactResidual(i, frame_position_[i]);
      const double distance = frame_position_[i].coeff(2);
      const double residual = data.slack.coeff(i) - dtau * distance;
      const double duality = computeDuality(data.slack.coeff(i), 
                                            data.dual.coeff(i));
      norm += residual * residual + duality * duality;
    }
  }
  return norm;
}


inline int DistanceToContactSurface::dimc_impl() const {
  return dimc_;
}


inline double DistanceToContactSurface::maxSlackStepSize_impl(
    const ConstraintComponentData& data,
    const std::vector<bool>& is_contact_active) const {
  double min_step_size = 1;
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      const double fraction_to_boundary 
          = fractionToBoundary(data.slack.coeff(i), data.dslack.coeff(i));
      if (fraction_to_boundary > 0 && fraction_to_boundary < 1) {
        if (fraction_to_boundary < min_step_size) {
          min_step_size = fraction_to_boundary;
        }
      }
    }
  }
  assert(min_step_size > 0);
  assert(min_step_size <= 1);
  return min_step_size;
}


inline double DistanceToContactSurface::maxDualStepSize_impl(
    const ConstraintComponentData& data,
    const std::vector<bool>& is_contact_active) const {
  double min_step_size = 1;
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      const double fraction_to_boundary 
          = fractionToBoundary(data.dual.coeff(i), data.ddual.coeff(i));
      if (fraction_to_boundary > 0 && fraction_to_boundary < 1) {
        if (fraction_to_boundary < min_step_size) {
          min_step_size = fraction_to_boundary;
        }
      }
    }
  }
  assert(min_step_size > 0);
  assert(min_step_size <= 1);
  return min_step_size;
}


inline void DistanceToContactSurface::updateSlack_impl(
    ConstraintComponentData& data, const std::vector<bool>& is_contact_active,
    const double step_size) const {
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      data.slack.coeffRef(i) += step_size * data.dslack.coeff(i);
    }
  }
}


inline void DistanceToContactSurface::updateDual_impl(
    ConstraintComponentData& data, const std::vector<bool>& is_contact_active,
    const double step_size) const {
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      data.dual.coeffRef(i) += step_size * data.ddual.coeff(i);
    }
  }
}


inline double DistanceToContactSurface::costSlackBarrier_impl(
    const ConstraintComponentData& data,
    const std::vector<bool>& is_contact_active) const {
  double cost = 0;
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      cost += costSlackBarrier(data.slack.coeff(i));
    }
  }
  return cost;
}


inline double DistanceToContactSurface::costSlackBarrier_impl(
    const ConstraintComponentData& data, 
    const std::vector<bool>& is_contact_active, const double step_size) const {
  double cost = 0;
  for (int i=0; i<dimc_; ++i) {
    if (!is_contact_active[i]) {
      cost += costSlackBarrier(data.slack.coeff(i), data.dslack.coeff(i), 
                               step_size);
    }
  }
  return cost;
}

} // namespace idocp

#endif // IDOCP_DISTANCE_TO_CONTACT_SURFACE_HXX_ 