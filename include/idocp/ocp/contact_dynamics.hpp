#ifndef IDOCP_CONTACT_DYNAMICS_HPP_
#define IDOCP_CONTACT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/contact_dynamics_data.hpp"


namespace idocp {

class ContactDynamics {
public:
  ContactDynamics(const Robot& robot);

  ContactDynamics();

  ~ContactDynamics();

  ContactDynamics(const ContactDynamics&) = default;

  ContactDynamics& operator=(const ContactDynamics&) = default;
 
  ContactDynamics(ContactDynamics&&) noexcept = default;

  ContactDynamics& operator=(ContactDynamics&&) noexcept = default;

  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const double dtau, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual);

  static void linearizeInverseDynamics(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const SplitSolution& s, 
                                       ContactDynamicsData& data);

  static void linearizeContactConstraint(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double dtau, 
                                         ContactDynamicsData& data);

  void condenseContactDynamics(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const double dtau, SplitKKTMatrix& kkt_matrix, 
                               SplitKKTResidual& kkt_residual);

  static void condensing(const Robot& robot, const double dtau, 
                         ContactDynamicsData& data, SplitKKTMatrix& kkt_matrix, 
                         SplitKKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const Robot& robot, SplitDirection& d);

  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, const double dtau, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     const SplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     SplitDirection& d);

  static void expansionPrimal(const Robot& robot, 
                              const ContactDynamicsData& data, 
                              SplitDirection& d);

  template <typename VectorType>
  static void expansionDual(const Robot& robot, const double dtau, 
                            ContactDynamicsData& data,
                            const SplitKKTMatrix& kkt_matrix, 
                            const SplitKKTResidual& kkt_residual,
                            const Eigen::MatrixBase<VectorType>& dgmm,
                            SplitDirection& d);

  void computeContactDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const double dtau, 
                                      const SplitSolution& s);

  double l1NormContactDynamicsResidual(const double dtau) const;

  double squaredNormContactDynamicsResidual(const double dtau) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ContactDynamicsData data_;
  bool has_floating_base_, has_active_contacts_;
  static constexpr int kDimFloatingBase = 6;

  void setContactStatus(const ContactStatus& contact_status);

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 