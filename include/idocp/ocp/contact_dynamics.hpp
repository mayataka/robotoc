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
#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {

///
/// @class ContactDynamics
/// @brief Contact dynamics constraint.
///
class ContactDynamics {
public:
  ///
  /// @brief Constructs the contact dynamics.
  /// @param[in] robot Robot model. 
  ///
  ContactDynamics(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ContactDynamics();

  ///
  /// @brief Destructor. 
  ///
  ~ContactDynamics();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactDynamics(const ContactDynamics&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ContactDynamics& operator=(const ContactDynamics&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactDynamics(ContactDynamics&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactDynamics& operator=(ContactDynamics&&) noexcept = default;

  ///
  /// @brief Linearizes the contact dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const double dt, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearizes the inverse dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] data Data for contact dynamics.
  ///
  static void linearizeInverseDynamics(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const SplitSolution& s, 
                                       ContactDynamicsData& data);

  ///
  /// @brief Linearizes the contact constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in, out] data Data for contact dynamics.
  ///
  static void linearizeContactConstraint(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         ContactDynamicsData& data);

  ///
  /// @brief Condenses the contact dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] is_forward_euler If true, the forward Euler is used. If false,
  /// the backward Euler is used.
  ///
  void condenseContactDynamics(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const double dt, SplitKKTMatrix& kkt_matrix, 
                               SplitKKTResidual& kkt_residual,
                               const bool is_forward_euler);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedPrimalDirection(const Robot& robot, 
                                       SplitDirection& d) const;

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in] dgmm Direction of the costate of the next time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, const double dt, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     const SplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     SplitDirection& d);

  ///
  /// @brief Condenses the switching constraint. 
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] jac Jacobian of the switching constraint.
  ///
  void condenseSwitchingConstraint(SplitKKTResidual& kkt_residual, 
                                   SplitStateConstraintJacobian& jac) const;

  ///
  /// @brief Computes the residual in the contact dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  ///
  void computeContactDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const SplitSolution& s);

  ///
  /// @brief Returns l1-norm of the residual in the contact dynamics constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return l1-norm of the residual in the contact dynamics constraint.
  ///
  double l1NormContactDynamicsResidual(const double dt) const;

  ///
  /// @brief Returns squared norm of the residual in the contact dynamics 
  /// constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return Squared norm of the residual in the contact dynamics constraint.
  ///
  double squaredNormContactDynamicsResidual(const double dt) const;

private:
  ContactDynamicsData data_;
  bool has_floating_base_, has_active_contacts_;
  int dimv_, dimu_, dim_passive_;
  static constexpr int kDimFloatingBase = 6;

  void setContactStatus(const ContactStatus& contact_status);

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 