#ifndef IDOCP_PDIPM_FUNCS_HPP_
#define IDOCP_PDIPM_FUNCS_HPP_

#include "Eigen/Core"


namespace idocp {
namespace pdipm {
namespace pdipmfunc {

void SetSlackAndDualPositive(const unsigned int dim, const double barrier,
                             Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual);

double FractionToBoundary(const unsigned int dim, const Eigen::VectorXd& vec, 
                          const Eigen::VectorXd& dvec);

void ComputeDualDirection(const double barrier, const Eigen::VectorXd& dual, 
                          const Eigen::VectorXd& slack, 
                          const Eigen::VectorXd& slack_direction, 
                          Eigen::VectorXd& dual_direction);

} // namespace pdipmfunc
} // namespace pdipm
} // namespace idocp


#endif // IDOCP_PDIPM_FUNCS_HPP_