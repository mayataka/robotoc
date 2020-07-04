#ifndef IDOCP_CONSTRAINTS_FB_FB_FUNC_HPP_
#define IDOCP_CONSTRAINTS_FB_FB_FUNC_HPP_

#include "Eigen/Core"


namespace idocp {
namespace fb {
namespace fbfunc {

void SetSlackAndDualPositive(const unsigned int dim, const double barrier,
                             Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual);

void ComputeFischerBurmeisterRadius(const unsigned int dim,
                                    const double barrier, 
                                    const Eigen::VectorXd& slack, 
                                    const Eigen::VectorXd& dual, 
                                    Eigen::VectorXd& radius);

double FractionToBoundary(const unsigned int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec);

void ComputeDualDirection(const Eigen::VectorXd& dual_tilde, 
                          const Eigen::VectorXd& slack_tilde, 
                          const Eigen::VectorXd& slack_direction, 
                          const Eigen::VectorXd& fb_residual, 
                          Eigen::VectorXd& dual_direction);

double SlackBarrierCost(const unsigned int dim, const double barrier, 
                        const Eigen::VectorXd& slack);

} // namespace fbfunc
} // namespace fb
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_FB_FB_FUNC_HPP_