#ifndef IDOCP_EQUALITY_CONSTRAINTS_HXX_
#define IDOCP_EQUALITY_CONSTRAINTS_HXX_

#include <assert.h>

namespace idocp {
namespace eqconstraints {

inline void AugmentEqualityConstraints(Robot& robot, const double dtau,
                                       const SplitSolution& s, 
                                       KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual) {
  assert(dtau > 0);
  if (robot.has_active_contacts()) {
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
    robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                      kkt_matrix.Ca());
    kkt_residual.lq() 
        += kkt_matrix.Cq().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
    kkt_residual.lv() 
        += kkt_matrix.Cv().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
    kkt_residual.la() 
        += kkt_matrix.Ca().topRows(robot.dimf()).transpose() 
            * s.mu.head(robot.dimf());
  }
  if (robot.has_floating_base()) {
    kkt_residual.C().tail(robot.dim_passive()) 
        = dtau * s.u.head(robot.dim_passive());
    kkt_residual.lu.head(robot.dim_passive()) 
        += dtau * s.mu_active().tail(robot.dim_passive());
  }
}


inline void AugmentCondensedEqualityConstraints(Robot& robot, const double dtau,
                                                const SplitSolution& s, 
                                                KKTMatrix& kkt_matrix, 
                                                KKTResidual& kkt_residual) {
  assert(dtau > 0);
  if (robot.has_active_contacts()) {
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
    robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                      kkt_matrix.Ca());
  }
  if (robot.has_active_contacts() || robot.has_floating_base()) {
    kkt_residual.lq().noalias() += kkt_matrix.Cq().transpose() * s.mu_active();
    kkt_residual.lv().noalias() += kkt_matrix.Cv().transpose() * s.mu_active();
    kkt_residual.la().noalias() += kkt_matrix.Ca().transpose() * s.mu_active();
  }
}


inline double ViolationL1Norm(const KKTResidual& kkt_residual) {
  return kkt_residual.C().lpNorm<1>();
}


inline double ComputeViolationL1Norm(Robot& robot, const double dtau, 
                                     const SplitSolution& s, 
                                     KKTResidual& kkt_residual) {
  assert(dtau > 0);
  if (robot.has_active_contacts()) {
    robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  }
  if (robot.has_floating_base()) {
    kkt_residual.C().tail(robot.dim_passive()) 
        = dtau * s.u.head(robot.dim_passive());
  }
  return kkt_residual.C().lpNorm<1>();
}

} // namespace eqconstraints 
} // namespace idocp 

#endif // IDOCP_EQUALITY_CONSTRAINTS_HXX_