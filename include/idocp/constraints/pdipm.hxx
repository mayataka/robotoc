#ifndef IDOCP_CONSTRAINTS_PDIPM_HXX_
#define IDOCP_CONSTRAINTS_PDIPM_HXX_

#include "idocp/constraints/pdipm.hpp"

#include <cmath>
#include <cassert>


namespace idocp {
namespace pdipm {

inline void SetSlackAndDualPositive(const double barrier, 
                                    ConstraintComponentData& data) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  for (int i=0; i<data.slack.size(); ++i) {
    while (data.slack.coeff(i) < barrier) {
      data.slack.coeffRef(i) += barrier;
    }
  }
  data.dual.array() = barrier / data.slack.array();
}


inline void ComputeComplementarySlackness(const double barrier, 
                                          ConstraintComponentData& data) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  data.cmpl.array() = data.slack.array() * data.dual.array() - barrier;
}


inline void ComputeComplementarySlackness(const double barrier, 
                                          ConstraintComponentData& data,
                                          const int start, const int size) {
  assert(barrier > 0);
  assert(data.checkDimensionalConsistency());
  data.cmpl.segment(start, size).array() 
      = data.slack.segment(start, size).array() 
          * data.dual.segment(start, size).array() - barrier;
}


inline double ComputeComplementarySlackness(const double barrier, 
                                            const double slack, 
                                            const double dual) {
  assert(barrier > 0);
  return (slack * dual - barrier); 
}


inline double FractionToBoundarySlack(const double fraction_rate, 
                                      const ConstraintComponentData& data) {
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(data.checkDimensionalConsistency());
  return FractionToBoundary(data.dimc(), fraction_rate, data.slack, data.dslack);
}


inline double FractionToBoundaryDual(const double fraction_rate, 
                                     const ConstraintComponentData& data) {
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(data.checkDimensionalConsistency());
  return FractionToBoundary(data.dimc(), fraction_rate, data.dual, data.ddual);
}


inline double FractionToBoundary(const int dim, const double fraction_rate, 
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


inline void ComputeDualDirection(ConstraintComponentData& data) {
  assert(data.checkDimensionalConsistency());
  data.ddual.array() 
      = - (data.dual.array()*data.dslack.array()+data.cmpl.array())
          / data.slack.array();
}


inline void ComputeDualDirection(ConstraintComponentData& data, 
                                 const int start, const int size) {
  data.ddual.segment(start, size).array() 
      = - (data.dual.segment(start, size).array()
              *data.dslack.segment(start, size).array()
            +data.cmpl.segment(start, size).array())
          / data.slack.segment(start, size).array();
}


inline double ComputeDualDirection(const double slack, const double dual,
                                   const double dslack, const double cmpl) {
  return (- (dual * dslack + cmpl) / slack);
}


template <typename VectorType>
inline double LogBarrier(const double barrier, 
                          const Eigen::MatrixBase<VectorType>& vec) {
  assert(barrier > 0);
  assert(vec.array().minCoeff() > 0);
  return (- barrier * vec.array().log().sum());
}

} // namespace pdipm
} // namespace idocp

#endif // IDOCP_CONSTRAINTS_PDIPM__HXX_