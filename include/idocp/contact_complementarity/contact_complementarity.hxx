#ifndef IDOCP_CONTACT_COMPLEMENTARITY_HXX_
#define IDOCP_CONTACT_COMPLEMENTARITY_HXX_

#include "idocp/contact_complementarity/contact_complementarity.hpp"

#include <algorithm>

namespace idocp {

inline ContactComplementarity::ContactComplementarity(const Robot& robot)
  : distance_to_contact_surface_constraint_(robot),
    contact_normal_force_constraint_(robot),
    friction_cone_constraint_(robot),
    distance_to_contact_surface_constraint_data_(
        distance_to_contact_surface_constraint_.dimc()),
    contact_normal_force_constraint_data_(contact_normal_force_constraint_.dimc()),
    friction_cone_constraint_data_(friction_cone_constraint_.dimc()),
    is_each_contact_active_(robot.max_point_contacts(), false) {
}


inline ContactComplementarity::ContactComplementarity()
  : distance_to_contact_surface_constraint_(),
    contact_normal_force_constraint_(),
    friction_cone_constraint_(),
    distance_to_contact_surface_constraint_data_(),
    contact_normal_force_constraint_data_(),
    friction_cone_constraint_data_(),
    is_each_contact_active_() {
}


inline ContactComplementarity::~ContactComplementarity() {
}


inline bool ContactComplementarity::isFeasible(Robot& robot, 
                                               const SplitSolution& s) {
  if (!is_each_contact_active_.empty()) {
    if (!distance_to_contact_surface_constraint_.isFeasible(
            robot, distance_to_contact_surface_constraint_data_, s)) {
      return false;
    }
    if (!contact_normal_force_constraint_.isFeasible(
            robot, contact_normal_force_constraint_data_, s)) {
      return false;
    }
    if (!friction_cone_constraint_.isFeasible(
            robot, friction_cone_constraint_data_, s)) {
      return false;
    }
  }
  return true;
}


inline void ContactComplementarity::setContactStatus(const Robot& robot) {
  is_each_contact_active_ = robot.is_each_contact_active();
}


inline void ContactComplementarity::setSlackAndDual(Robot& robot, 
                                                    const double dtau, 
                                                    const SplitSolution& s) {
  assert(dtau > 0);
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.setSlackAndDual(
        robot, distance_to_contact_surface_constraint_data_, dtau, s);
    contact_normal_force_constraint_.setSlackAndDual(
        robot, contact_normal_force_constraint_data_, dtau, s);
    friction_cone_constraint_.setSlackAndDual(
        robot, friction_cone_constraint_data_, dtau, s);
  }
}


inline void ContactComplementarity::augmentDualResidual(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTResidual& kkt_residual) {
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.augmentDualResidual(
        robot, distance_to_contact_surface_constraint_data_, dtau, s,
        kkt_residual);
    contact_normal_force_constraint_.augmentDualResidual(
        robot, contact_normal_force_constraint_data_, dtau, s, kkt_residual);
    friction_cone_constraint_.augmentDualResidual(
        robot, friction_cone_constraint_data_, dtau, s, kkt_residual);
  }
}


inline void ContactComplementarity::condenseSlackAndDual(
    Robot& robot, const double dtau, const SplitSolution& s, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.condenseSlackAndDual(
        robot, distance_to_contact_surface_constraint_data_, dtau, s, 
        kkt_matrix, kkt_residual);
    contact_normal_force_constraint_.condenseSlackAndDual(
        robot, contact_normal_force_constraint_data_, dtau, s, 
        kkt_matrix, kkt_residual);
    friction_cone_constraint_.condenseSlackAndDual(
        robot, friction_cone_constraint_data_, dtau, s, 
        kkt_matrix, kkt_residual);
  }
}


inline void ContactComplementarity::computeSlackAndDualDirection(
    Robot& robot, const double dtau, const SplitSolution& s, 
    const SplitDirection& d) {
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.computeSlackAndDualDirection(
        robot, distance_to_contact_surface_constraint_data_, dtau, s, d);
    contact_normal_force_constraint_.computeSlackAndDualDirection(
        robot, contact_normal_force_constraint_data_, dtau, s, d);
    friction_cone_constraint_.computeSlackAndDualDirection(
        robot, friction_cone_constraint_data_, dtau, s, d);
  }
}


inline double ContactComplementarity::residualL1Nrom(
    Robot& robot, const double dtau, const SplitSolution& s) {
  double norm = 0;
  if (!is_each_contact_active_.empty()) {
    norm += distance_to_contact_surface_constraint_.residualL1Nrom(
        robot, distance_to_contact_surface_constraint_data_, dtau, s);
    norm += contact_normal_force_constraint_.residualL1Nrom(
        robot, contact_normal_force_constraint_data_, dtau, s);
    norm += friction_cone_constraint_.residualL1Nrom(
        robot, friction_cone_constraint_data_, dtau, s);
  }
  return norm;
}


inline double ContactComplementarity::squaredKKTErrorNorm(
    Robot& robot, const double dtau, const SplitSolution& s) {
  double norm = 0;
  if (!is_each_contact_active_.empty()) {
    norm += distance_to_contact_surface_constraint_.squaredKKTErrorNorm(
        robot, distance_to_contact_surface_constraint_data_, dtau, s);
    norm += contact_normal_force_constraint_.squaredKKTErrorNorm(
        robot, contact_normal_force_constraint_data_, dtau, s);
    norm += friction_cone_constraint_.squaredKKTErrorNorm(
        robot, friction_cone_constraint_data_, dtau, s);
  }
  return norm;
}


inline double ContactComplementarity::maxSlackStepSize() const {
  double norm = 0;
  if (!is_each_contact_active_.empty()) {
    const double step_size1 
        = distance_to_contact_surface_constraint_.maxSlackStepSize(
              distance_to_contact_surface_constraint_data_,
              is_each_contact_active_);
    const double step_size2 
        = contact_normal_force_constraint_.maxSlackStepSize(
              contact_normal_force_constraint_data_, is_each_contact_active_);
    const double step_size3 
        = friction_cone_constraint_.maxSlackStepSize(
              friction_cone_constraint_data_, is_each_contact_active_);
    return std::min({step_size1, step_size2, step_size3});
  }
  return 1;
}


inline double ContactComplementarity::maxDualStepSize() const {
  double norm = 0;
  if (!is_each_contact_active_.empty()) {
    const double step_size1 
        = distance_to_contact_surface_constraint_.maxDualStepSize(
              distance_to_contact_surface_constraint_data_, 
              is_each_contact_active_);
    const double step_size2 
        = contact_normal_force_constraint_.maxDualStepSize(
              contact_normal_force_constraint_data_, is_each_contact_active_);
    const double step_size3 
        = friction_cone_constraint_.maxDualStepSize(
              friction_cone_constraint_data_, is_each_contact_active_);
    return std::min({step_size1, step_size2, step_size3});
  }
  return 1;
}


inline void ContactComplementarity::updateSlack(const double step_size) {
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.updateSlack(
        distance_to_contact_surface_constraint_data_, is_each_contact_active_, 
        step_size);
    contact_normal_force_constraint_.updateSlack(
        contact_normal_force_constraint_data_, is_each_contact_active_, 
        step_size);
    friction_cone_constraint_.updateSlack(
        friction_cone_constraint_data_, is_each_contact_active_, step_size);
  }
}


inline void ContactComplementarity::updateDual(const double step_size) {
  if (!is_each_contact_active_.empty()) {
    distance_to_contact_surface_constraint_.updateDual(
        distance_to_contact_surface_constraint_data_, is_each_contact_active_, 
        step_size);
    contact_normal_force_constraint_.updateDual(
        contact_normal_force_constraint_data_, is_each_contact_active_, 
        step_size);
    friction_cone_constraint_.updateDual(
        friction_cone_constraint_data_, is_each_contact_active_, step_size);
  }
}


inline double ContactComplementarity::costSlackBarrier() const {
  double barrier = 0;
  if (!is_each_contact_active_.empty()) {
    barrier += distance_to_contact_surface_constraint_.costSlackBarrier(
        distance_to_contact_surface_constraint_data_, is_each_contact_active_);
    barrier += contact_normal_force_constraint_.costSlackBarrier(
        contact_normal_force_constraint_data_, is_each_contact_active_);
    barrier += friction_cone_constraint_.costSlackBarrier(
        friction_cone_constraint_data_, is_each_contact_active_);
  }
  return barrier;
}


inline double ContactComplementarity::costSlackBarrier(
    const double step_size) const {
  double barrier = 0;
  if (!is_each_contact_active_.empty()) {
    barrier += distance_to_contact_surface_constraint_.costSlackBarrier(
        distance_to_contact_surface_constraint_data_, is_each_contact_active_, 
        step_size);
    barrier += contact_normal_force_constraint_.costSlackBarrier(
        contact_normal_force_constraint_data_, is_each_contact_active_, 
        step_size);
    barrier += friction_cone_constraint_.costSlackBarrier(
        friction_cone_constraint_data_, is_each_contact_active_, step_size);
  }
  return barrier;
}


inline void ContactComplementarity::setBarrier(const double barrier) {
  distance_to_contact_surface_constraint_.setBarrier(barrier);
  contact_normal_force_constraint_.setBarrier(barrier);
  friction_cone_constraint_.setBarrier(barrier);
}


inline void ContactComplementarity::setFractionToBoundaryRate(
    const double fraction_to_boundary_rate) {
  distance_to_contact_surface_constraint_.setFractionToBoundaryRate(
      fraction_to_boundary_rate);
  contact_normal_force_constraint_.setFractionToBoundaryRate(
      fraction_to_boundary_rate);
  friction_cone_constraint_.setFractionToBoundaryRate(
      fraction_to_boundary_rate);
}

} // namespace idocp

#endif // IDOCP_CONTACT_COMPLEMENTARITY_HXX_ 