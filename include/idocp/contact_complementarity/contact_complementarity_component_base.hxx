#ifndef IDOCP_CONTACT_COMPLEMENTARITY_COMPONENT_BASE_HXX_ 
#define IDOCP_CONTACT_COMPLEMENTARITY_COMPONENT_BASE_HXX_

#include "idocp/contact_complementarity/contact_complementarity_component_base.hpp"
#include "idocp/constraints/pdipm_func.hpp"

#include <cmath>
#include <exception>
#include <assert.h>

namespace idocp {

template <typename Derived>
inline ContactComplementarityComponentBase<Derived>::
ContactComplementarityComponentBase(const double barrier, 
                                    const double fraction_to_boundary_rate) 
  : barrier_(barrier),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
    try {
    if (barrier <= 0) {
      throw std::out_of_range(
          "invalid argment: barrirer must be positive");
    }
    if (fraction_to_boundary_rate <= 0) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rate must be positive");
    }
    if (fraction_to_boundary_rate >= 1) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rate must be less than 1");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


template <typename Derived>
inline ContactComplementarityComponentBase<Derived>::
ContactComplementarityComponentBase()
  : barrier_(0),
    fraction_to_boundary_rate_(0) {
}


template <typename Derived>
inline ContactComplementarityComponentBase<Derived>::
~ContactComplementarityComponentBase() {
}


template <typename Derived>
inline bool ContactComplementarityComponentBase<Derived>::isFeasible(
    Robot& robot, ConstraintComponentData& data, 
    const SplitSolution& s) const {
  return static_cast<const Derived*>(this)->isFeasible_impl(robot, data, s);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  static_cast<const Derived*>(this)->setSlackAndDual_impl(robot, data, dtau, s);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  static_cast<const Derived*>(this)->augmentDualResidual_impl(robot, data, dtau,  
                                                              s, kkt_residual);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  static_cast<const Derived*>(this)->condenseSlackAndDual_impl(robot, data, dtau, 
                                                               s, kkt_matrix, 
                                                               kkt_residual);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::
computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                             const double dtau, const SplitSolution& s, 
                             const SplitDirection& d) const {
  static_cast<const Derived*>(this)->computeSlackAndDualDirection_impl(robot, data, 
                                                                 dtau, s, d);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, 
    const double dtau, const SplitSolution& s) const {
  return static_cast<const Derived*>(this)->residualL1Nrom_impl(robot, data, dtau, s);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  return static_cast<const Derived*>(this)->squaredKKTErrorNorm_impl(robot, data, 
                                                               dtau, s);
}


template <typename Derived>
inline int ContactComplementarityComponentBase<Derived>::dimc() const {
  return static_cast<const Derived*>(this)->dimc_impl();
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::maxSlackStepSize(
    const ConstraintComponentData& data, 
    const std::vector<bool>& is_contact_active) const {
  return static_cast<const Derived*>(this)->maxSlackStepSize_impl(
      data, is_contact_active);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::maxDualStepSize(
    const ConstraintComponentData& data, 
    const std::vector<bool>& is_contact_active) const {
  return static_cast<const Derived*>(this)->maxDualStepSize_impl(
      data, is_contact_active);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::updateSlack(
    ConstraintComponentData& data, const std::vector<bool>& is_contact_active,
    const double step_size) const {
  assert(step_size > 0);
  static_cast<const Derived*>(this)->updateSlack_impl(data, is_contact_active, step_size);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::updateDual(
    ConstraintComponentData& data, const std::vector<bool>& is_contact_active,
    const double step_size) const {
  assert(step_size > 0);
  static_cast<const Derived*>(this)->updateDual_impl(data, is_contact_active, step_size);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::costSlackBarrier(
    const ConstraintComponentData& data, 
    const std::vector<bool>& is_contact_active) const {
  return static_cast<const Derived*>(this)->costSlackBarrier_impl(
      data, is_contact_active);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::costSlackBarrier(
    const ConstraintComponentData& data, 
    const std::vector<bool>& is_contact_active,
    const double step_size) const {
  return static_cast<const Derived*>(this)->costSlackBarrier_impl(
      data, is_contact_active, step_size);
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::setBarrier(
    const double barrier) {
  assert(barrier > 0);
  barrier_ = barrier;
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::
setFractionToBoundaryRate(const double fraction_to_boundary_rate) {
  assert(fraction_to_boundary_rate > 0);
  fraction_to_boundary_rate_ = fraction_to_boundary_rate;
}


template <typename Derived>
inline void ContactComplementarityComponentBase<Derived>::
setSlackAndDualPositive(Eigen::VectorXd& slack, Eigen::VectorXd& dual) const {
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack, dual);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::costSlackBarrier(
    const double slack) const {
  return - barrier_ * std::log(slack);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::costSlackBarrier(
    const double slack, const double dslack, const double step_size) const {
  return - barrier_ * std::log(slack + step_size * dslack);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::computeDuality(
    const double slack, const double dual) const {
  return pdipmfunc::ComputeDuality(barrier_, slack, dual);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double duality) const {
  return pdipmfunc::ComputeDualDirection(slack, dual, dslack, duality);
}


template <typename Derived>
inline double ContactComplementarityComponentBase<Derived>::fractionToBoundary(
    const double vec, const double dvec) const {
  return pdipmfunc::FractionToBoundary(fraction_to_boundary_rate_, vec, dvec);
}

} // namespace idocp

#endif // IDOCP_CONTACT_COMPLEMENTARITY_COMPONENT_BASE_HXX_ 