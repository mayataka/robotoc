#include "idocp/constraints/pdipm_func.hpp"

#include <assert.h>


namespace idocp {
namespace pdipmfunc {

inline void SetSlackAndDualPositive(const double barrier, 
                                    Eigen::VectorXd& slack, 
                                    Eigen::VectorXd& dual) {
  assert(barrier > 0);
  assert(slack.size() == dual.size());
  for (int i=0; i<slack.size(); ++i) {
    while (slack.coeff(i) < barrier) {
      slack.coeffRef(i) += barrier;
    }
  }
  dual.array() = barrier / slack.array();
}


inline void ComputeDuality(const double barrier, const Eigen::VectorXd& slack, 
                           const Eigen::VectorXd& dual, 
                           Eigen::VectorXd& duality) {
  assert(barrier > 0);
  assert(slack.size() == dual.size());
  assert(slack.size() == duality.size());
  duality.array() = slack.array() * dual.array() - barrier;
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
    double fraction_to_boundary = - fraction_rate * (vec.coeff(i)/dvec.coeff(i));
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


inline void ComputeDualDirection(const Eigen::VectorXd& slack, 
                                 const Eigen::VectorXd& dual, 
                                 const Eigen::VectorXd& slack_direction, 
                                 const Eigen::VectorXd& duality, 
                                 Eigen::VectorXd& dual_direction) {
  assert(slack.size() == dual.size());
  assert(slack.size() == slack_direction.size());
  assert(slack.size() == duality.size());
  assert(slack.size() == dual_direction.size());
  dual_direction.array() 
      = - (dual.array()*slack_direction.array()+duality.array())
          / slack.array();
}


inline double CostSlackBarrier(const double barrier, 
                               const Eigen::VectorXd& slack) {
  assert(barrier > 0);
  const double cost = - barrier * slack.array().log().sum();
  return cost;
}

} // namespace pdipm
} // namespace idocp
