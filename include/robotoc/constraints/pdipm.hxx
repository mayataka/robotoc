#ifndef ROBOTOC_CONSTRAINTS_PDIPM_HXX_
#define ROBOTOC_CONSTRAINTS_PDIPM_HXX_

#include "robotoc/constraints/pdipm.hpp"

#include <cmath>
#include <cassert>


namespace robotoc {
namespace pdipm {

inline void setSlackAndDualPositive(const double barrier, 
                                    ConstraintComponentData& data) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  const double sqrt_barrier = std::sqrt(barrier);
  for (int i=0; i<data.slack.size(); ++i) {
    if (data.slack.coeff(i) < sqrt_barrier) {
      data.slack.coeffRef(i) = sqrt_barrier;
    }
  }
  data.dual.array() = barrier / data.slack.array();
}


inline void computeComplementarySlackness(const double barrier, 
                                          ConstraintComponentData& data) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  data.cmpl.array() = data.slack.array() * data.dual.array() - barrier;
}


inline void computeComplementarySlackness(const double barrier, 
                                          ConstraintComponentData& data,
                                          const int start, const int size) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  data.cmpl.segment(start, size).array() 
      = data.slack.segment(start, size).array() 
          * data.dual.segment(start, size).array() - barrier;
}


template <int Size>
inline void computeComplementarySlackness(const double barrier, 
                                          ConstraintComponentData& data,
                                          const int start) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  data.cmpl.template segment<Size>(start).array() 
      = data.slack.template segment<Size>(start).array() 
          * data.dual.template segment<Size>(start).array() - barrier;
}


inline double computeComplementarySlackness(const double barrier, 
                                            const double slack, 
                                            const double dual) {
  assert(barrier > 0);
  return (slack * dual - barrier); 
}


inline void computeCondensingCoeffcient(ConstraintComponentData& data) {
  assert(data.checkDimensionalConsistency());
  data.cond.array() = (data.dual.array()*data.residual.array()-data.cmpl.array()) 
                        / data.slack.array();
}


inline void computeCondensingCoeffcient(ConstraintComponentData& data,
                                        const int start, const int size) {
  assert(data.checkDimensionalConsistency());
  data.cond.segment(start, size).array() 
      = (data.dual.segment(start, size).array()
          *data.residual.segment(start, size).array()
          -data.cmpl.segment(start, size).array()) 
         / data.slack.segment(start, size).array();
}


template <int Size>
inline void computeCondensingCoeffcient(ConstraintComponentData& data,
                                        const int start) {
  assert(data.checkDimensionalConsistency());
  data.cond.template segment<Size>(start).array() 
      = (data.dual.template segment<Size>(start).array()
          *data.residual.template segment<Size>(start).array()
          -data.cmpl.template segment<Size>(start).array()) 
         / data.slack.template segment<Size>(start).array();
}


inline double computeCondensingCoeffcient(const double slack, const double dual,
                                          const double residual, 
                                          const double cmpl) {
  return ((dual*residual-cmpl)/slack);
}


inline double fractionToBoundarySlack(const double fraction_rate, 
                                      const ConstraintComponentData& data) {
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(data.checkDimensionalConsistency());
  return fractionToBoundary(data.dimc(), fraction_rate, data.slack, data.dslack);
}


inline double fractionToBoundaryDual(const double fraction_rate, 
                                     const ConstraintComponentData& data) {
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(data.checkDimensionalConsistency());
  return fractionToBoundary(data.dimc(), fraction_rate, data.dual, data.ddual);
}


inline double fractionToBoundary(const int dim, const double fraction_rate, 
                                 const Eigen::VectorXd& vec, 
                                 const Eigen::VectorXd& dvec) {
  assert(dim > 0);
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(vec.size() == dim);
  assert(dvec.size() == dim);
  double min_fraction_to_boundary = 1;
  for (int i=0; i<dim; ++i) {
    const double fraction_to_boundary 
        = - fraction_rate * (vec.coeff(i)/dvec.coeff(i));
    if (fraction_to_boundary > 0 && fraction_to_boundary < 1) {
      if (fraction_to_boundary < min_fraction_to_boundary) {
        min_fraction_to_boundary = fraction_to_boundary;
      }
    }
  }
  assert(min_fraction_to_boundary > 0);
  assert(min_fraction_to_boundary <= 1);
  return min_fraction_to_boundary;
}


inline double fractionToBoundary(const double fraction_rate, 
                                 const double var, const double dvar) {
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  const double fraction_to_boundary = - fraction_rate * (var/dvar);
  if (fraction_to_boundary > 0 && fraction_to_boundary < 1) {
    return fraction_to_boundary;
  }
  else {
    return 1.0;
  }
}


inline void computeDualDirection(ConstraintComponentData& data) {
  assert(data.checkDimensionalConsistency());
  data.ddual.array() 
      = - (data.dual.array()*data.dslack.array()+data.cmpl.array())
          / data.slack.array();
}


inline void computeDualDirection(ConstraintComponentData& data, 
                                 const int start, const int size) {
  data.ddual.segment(start, size).array() 
      = - (data.dual.segment(start, size).array()
              *data.dslack.segment(start, size).array()
            +data.cmpl.segment(start, size).array())
          / data.slack.segment(start, size).array();
}


template <int Size>
inline void computeDualDirection(ConstraintComponentData& data, 
                                 const int start) {
  data.ddual.template segment<Size>(start).array() 
      = - (data.dual.template segment<Size>(start).array()
              *data.dslack.template segment<Size>(start).array()
            +data.cmpl.template segment<Size>(start).array())
          / data.slack.template segment<Size>(start).array();
}


inline double computeDualDirection(const double slack, const double dual,
                                   const double dslack, const double cmpl) {
  return (- (dual * dslack + cmpl) / slack);
}


template <typename VectorType>
inline double logBarrier(const double barrier, 
                          const Eigen::MatrixBase<VectorType>& vec) {
  assert(barrier > 0);
  assert(vec.array().minCoeff() > 0);
  return (- barrier * vec.array().log().sum());
}

} // namespace pdipm
} // namespace robotoc

#endif // ROBOTOC_CONSTRAINTS_PDIPM__HXX_