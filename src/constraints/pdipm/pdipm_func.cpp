#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {
namespace pdipmfunc {


void SetSlackAndDualPositive(const unsigned int dim, const double barrier,
                             Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual) {
  assert(dim > 0);
  assert(slack.size() == dim);
  assert(dual.size() == dim);
  for (int i=0; i<dim; ++i) {
    while (slack.coeff(i) < barrier) {
      slack.coeffRef(i) += barrier;
    }
  }
  dual.array() = barrier / slack.array();
}


double FractionToBoundary(const unsigned int dim, const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec) {
  assert(dim > 0);
  assert(vec.size() == dim);
  assert(dvec.size() == dim);
  double max_fraction_to_boundary = 1;
  for (int i=0; i<dim; ++i) {
    const double fraction_to_boundary = (vec.coeff(i)/dvec.coeff(i));
    if (fraction_to_boundary > 0) {
      if (fraction_to_boundary < max_fraction_to_boundary) {
        max_fraction_to_boundary = fraction_to_boundary;
      }
    }
  }
  assert(max_fraction_to_boundary <= 1);
  assert(max_fraction_to_boundary > 0);
  return max_fraction_to_boundary;
}


void ComputeDualDirection(const double barrier, const Eigen::VectorXd& dual, 
                          const Eigen::VectorXd& slack, 
                          const Eigen::VectorXd& slack_direction, 
                          Eigen::VectorXd& dual_direction) {
  assert(barrier > 0);
  assert(dual.size() == slack.size());
  assert(dual.size() == slack_direction.size());
  assert(dual.size() == dual_direction.size());
  dual_direction.array() 
      = (-dual.array()*(slack.array()+slack_direction.array())+barrier)
        / slack.array();
}


} // namespace pdipmfunc
} // namespace pdipm
} // namespace idocp