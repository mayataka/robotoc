#ifndef IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_
#define IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_

#include "Eigen/Core"


namespace idocp {
namespace pdipm {
namespace pdipmfunc {

void SetSlackAndDualPositive(const int dim, const double barrier,
                             Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual);

void ComputeDualityResidual(const double barrier, const Eigen::VectorXd& slack, 
                            const Eigen::VectorXd& dual, 
                            Eigen::VectorXd& duality_residual);

double FractionToBoundary(const int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec);

void ComputeDualDirection(const Eigen::VectorXd& dual, 
                          const Eigen::VectorXd& slack, 
                          const Eigen::VectorXd& slack_direction, 
                          const Eigen::VectorXd& duality, 
                          Eigen::VectorXd& dual_direction);

double SlackBarrierCost(const int dim, const double barrier, 
                        const Eigen::VectorXd& slack);

} // namespace pdipmfunc
} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_