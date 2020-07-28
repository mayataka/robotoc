#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {
namespace pdipmfunc {


void SetSlackAndDualPositive(const int dim, const double barrier,
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


void ComputeDualityResidual(const double barrier, const Eigen::VectorXd& slack, 
                            const Eigen::VectorXd& dual, 
                            Eigen::VectorXd& duality_residual) {
  assert(barrier > 0);
  assert(slack.size() == dual.size());
  assert(slack.size() == duality_residual.size());
  duality_residual.array() = slack.array() * dual.array() - barrier;
}


double FractionToBoundary(const int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec) {
  assert(dim > 0);
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(vec.size() == dim);
  assert(dvec.size() == dim);
  double min_fraction_to_boundary = 1;
  for (int i=0; i<dim; ++i) {
    double fraction_to_boundary = - 0.995 * (vec.coeff(i)/dvec.coeff(i));
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


void ComputeDualDirection(const Eigen::VectorXd& dual, 
                          const Eigen::VectorXd& slack, 
                          const Eigen::VectorXd& slack_direction, 
                          const Eigen::VectorXd& duality, 
                          Eigen::VectorXd& dual_direction) {
  assert(dual.size() == slack.size());
  assert(dual.size() == slack_direction.size());
  assert(dual.size() == duality.size());
  assert(dual.size() == dual_direction.size());
  dual_direction.array() 
      = - (dual.array()*slack_direction.array()+duality.array())
          / slack.array();
}


double SlackBarrierCost(const int dim, const double barrier, 
                        const Eigen::VectorXd& slack) {
  assert(dim > 0);
  assert(barrier > 0);
  assert(slack.size() == dim);
  double cost = - barrier * slack.array().log().sum();
  return cost;
}

} // namespace pdipmfunc
} // namespace pdipm
} // namespace idocp