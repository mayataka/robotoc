#ifndef IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_
#define IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_

#include "Eigen/Core"


namespace idocp {
namespace pdipmfunc {

void SetSlackAndDualPositive(const double barrier, Eigen::VectorXd& slack, 
                             Eigen::VectorXd& dual);

void ComputeDuality(const double barrier, const Eigen::VectorXd& slack, 
                    const Eigen::VectorXd& dual, Eigen::VectorXd& duality);

double FractionToBoundary(const int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec,
                          const Eigen::VectorXd& dvec);

void ComputeDualDirection(const Eigen::VectorXd& slack, 
                          const Eigen::VectorXd& dual, 
                          const Eigen::VectorXd& slack_direction, 
                          const Eigen::VectorXd& duality, 
                          Eigen::VectorXd& dual_direction);

double CostSlackBarrier(const double barrier, const Eigen::VectorXd& slack);

double ComputeDuality(const double barrier, const double slack, const double dual);

double FractionToBoundary(const double fraction_rate, const double var,
                          const double dvar);

double ComputeDualDirection(const double slack, const double dual, 
                            const double slack_direction, const double duality);

double CostSlackBarrier(const double barrier, const double slack);

} // namespace pdipm
} // namespace idocp

#include "idocp/constraints/pdipm_func.hxx"

#endif // IDOCP_CONSTRAINTS_PDIPM_FUNC_HPP_