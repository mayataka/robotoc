#ifndef IDOCP_EQUALITY_CONSTRAINTS_HPP_
#define IDOCP_EQUALITY_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {
namespace eqconstraints {

void AugmentEqualityConstraints(Robot& robot, const double dtau, 
                                const SplitSolution& s, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual);

void AugmentCondensedEqualityConstraints(Robot& robot, const double dtau,
                                         const SplitSolution& s, 
                                         KKTMatrix& kkt_matrix, 
                                         KKTResidual& kkt_residual);

double ViolationL1Norm(const KKTResidual& kkt_residual);

double ComputeViolationL1Norm(Robot& robot, const double dtau, 
                              const SplitSolution& s, 
                              KKTResidual& kkt_residual);

} // namespace eqconstraints 
} // namespace idocp 

#include "idocp/ocp/equality_constraints.hxx"

#endif // IDOCP_EQUALITY_CONSTRAINTS_HPP_