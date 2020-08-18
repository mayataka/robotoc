#ifndef IDOCP_EQUALITY_CONSTRAINTS_HXX_
#define IDOCP_EQUALITY_CONSTRAINTS_HXX_

namespace idocp {
namespace equalityconstraints {

inline void LinearizeEqualityConstraints(Robot& robot, const double dtau,
                                         const SplitSolution& s, 
                                         KKTMatrix& kkt_matrix, 
                                         KKTResidual& kkt_residual) {
  LinearizeContactConstraints(robot, dtau, s, kkt_matrix, kkt_residual);
  LinearizeFloatingBaseConstraints(robot, dtau, s, kkt_residual);
}


inline void LinearizeContactConstraints(Robot& robot, const double dtau,
                                        const SplitSolution& s, 
                                        KKTMatrix& kkt_matrix, 
                                        KKTResidual& kkt_residual) {
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
}


inline void LinearizeFloatingBaseConstraints(const Robot& robot, 
                                             const double dtau,
                                             const SplitSolution& s, 
                                             KKTResidual& kkt_residual) {
  if (robot.has_floating_base()) {
    kkt_residual.C().tail(robot.dim_passive()) 
        = dtau * s.u.head(robot.dim_passive());
    kkt_residual.lu.head(robot.dim_passive())
        += dtau * s.mu_active().tail(robot.dim_passive());
  }
}


inline double ViolationL1Norm(const KKTResidual& kkt_residual) {
  return kkt_residual.C().lpNorm<1>();
}


inline double ViolationL1Norm(Robot& robot, const double dtau, 
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

} // namespace equalityconstraints 
} // namespace idocp 

#endif // IDOCP_EQUALITY_CONSTRAINTS_HXX_