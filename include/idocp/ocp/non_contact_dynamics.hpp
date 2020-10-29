#ifndef IDOCP_NON_CONTACT_DYNAMICS_HPP_
#define IDOCP_NON_CONTACT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class NonContactDynamics {
public:
  NonContactDynamics(const Robot& robot);

  NonContactDynamics();

  ~NonContactDynamics();

  NonContactDynamics(const NonContactDynamics&) = default;

  NonContactDynamics& operator=(const NonContactDynamics&) = default;
 
  NonContactDynamics(NonContactDynamics&&) noexcept = default;

  NonContactDynamics& operator=(NonContactDynamics&&) noexcept = default;

  void linearizeNonContactDynamics(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   const double dtau, const SplitSolution& s, 
                                   KKTMatrix& kkt_matrix, 
                                   KKTResidual& kkt_residual);

  void condenseNonContactDynamics(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const double dtau, const SplitSolution& s, 
                                  KKTMatrix& kkt_matrix, 
                                  KKTResidual& kkt_residual);
  
  void computeCondensedDirection(const double dtau, 
                                 const KKTMatrix& kkt_matrix, 
                                 const KKTResidual& kkt_residual, 
                                 SplitDirection& d);

  void computeNonContactDynamicsResidual(Robot& robot, 
                                         const ContactStatus& contact_status,
                                         const double dtau, 
                                         const SplitSolution& s, 
                                         KKTResidual& kkt_residual);

  double l1NormNonContactDynamicsResidual(const double dtau) const;

  double squaredNormNonContactDynamicsResidual(const double dtau) const;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kaq,
                            const Eigen::MatrixBase<MatrixType2>& Kav,
                            const Eigen::MatrixBase<MatrixType3>& Kuq,
                            const Eigen::MatrixBase<MatrixType4>& Kuv) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd ID_, lu_condensed_;
  Eigen::MatrixXd dID_dq_, dID_dv_, dID_da_, Quu_, Quu_du_dq_, Quu_du_dv_, 
                  Quu_du_da_;
  int dimv_;

  void linearizeInverseDynamics(Robot& robot, 
                                const ContactStatus& contact_status,
                                const SplitSolution& s, 
                                KKTResidual& kkt_residual);

  static void setContactForcesZero(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   const SplitSolution& s);

  static void computeInverseDynamicsResidual(Robot& robot, 
                                             const SplitSolution& s, 
                                             KKTResidual& kkt_residual);

};

} // namespace idocp 

#include "idocp/ocp/non_contact_dynamics.hxx"

#endif // IDOCP_NON_CONTACT_DYNAMICS_HPP_ 