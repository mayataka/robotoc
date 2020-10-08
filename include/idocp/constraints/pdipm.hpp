#ifndef IDOCP_CONSTRAINTS_PDIPM_HPP_
#define IDOCP_CONSTRAINTS_PDIPM_HPP_

#include "Eigen/Core"

#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {
namespace pdipm {

void SetSlackAndDualPositive(const double barrier, 
                             ConstraintComponentData& data);

void ComputeDuality(const double barrier, ConstraintComponentData& data);

double FractionToBoundarySlack(const double fraction_rate, 
                               const ConstraintComponentData& data);

double FractionToBoundaryDual(const double fraction_rate, 
                              const ConstraintComponentData& data);

double FractionToBoundary(const int dim, const double fraction_rate, 
                          const Eigen::VectorXd& vec,
                          const Eigen::VectorXd& dvec);

void ComputeDualDirection(ConstraintComponentData& data);

} // namespace pdipm
} // namespace idocp

#include "idocp/constraints/pdipm.hxx"

#endif // IDOCP_CONSTRAINTS_PDIPM_HPP_