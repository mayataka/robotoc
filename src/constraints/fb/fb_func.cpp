#include "constraints/fb/fb_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace fb {
namespace fbfunc {


void SetSlackAndDualPositive(const unsigned int dim, const double barrier,
                             Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual) {
  assert(dim > 0);
  assert(slack.size() == dim);
  assert(dual.size() == dim);
  for (unsigned int i=0; i<dim; ++i) {
    while (slack.coeff(i) < barrier) {
      slack.coeffRef(i) += barrier;
    }
  }
  dual.array() = barrier / slack.array();
}


void ComputeFischerBurmeisterRadius(const unsigned int dim, 
                                    const double barrier, 
                                    const Eigen::VectorXd& slack, 
                                    const Eigen::VectorXd& dual, 
                                    Eigen::VectorXd& radius) {
  assert(dim > 0);
  assert(barrier > 0);
  assert(slack.size() == dim);
  assert(dual.size() == dim);
  assert(radius.size() == dim);
  for (unsigned int i=0; i<dim; ++i) {
    radius.coeffRef(i) = std::sqrt(slack.coeff(i)*slack.coeff(i)
                                   +dual.coeff(i)*dual.coeff(i) 
                                   +2*barrier*barrier);
  }
}


double FractionToBoundary(const unsigned int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec) {
  assert(dim > 0);
  assert(fraction_rate > 0);
  assert(fraction_rate <= 1);
  assert(vec.size() == dim);
  assert(dvec.size() == dim);
  double min_fraction_to_boundary = 1;
  for (unsigned int i=0; i<dim; ++i) {
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


void ComputeDualDirection(const Eigen::VectorXd& dual_tilde, 
                          const Eigen::VectorXd& slack_tilde, 
                          const Eigen::VectorXd& slack_direction, 
                          const Eigen::VectorXd& fb_residual, 
                          Eigen::VectorXd& dual_direction) {
  dual_direction.array() 
      = (-slack_tilde.array()*slack_direction.array()+fb_residual.array())
          / dual_tilde.array();
}


double SlackBarrierCost(const unsigned int dim, const double barrier, 
                        const Eigen::VectorXd& slack) {
  assert(dim > 0);
  assert(barrier > 0);
  assert(slack.size() == dim);
  double cost = - barrier * slack.array().log().sum();
  return cost;
}

} // namespace fbfunc
} // namespace fb
} // namespace idocp