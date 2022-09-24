#ifndef ROBOTOC_CONTACT_DYNAMICS_HPP_
#define ROBOTOC_CONTACT_DYNAMICS_HPP_

#include <limits>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"


namespace robotoc {

///
/// @brief Computes the residual in the contact dynamics constraint. 
/// @param[in] robot Robot model. 
/// @param[in] contact_status Contact status of this time stage. 
/// @param[in] s Split solution of this time stage.
/// @param[in, out] data Data structure for the contact dynamics.
///
void evalContactDynamics(Robot& robot, const ContactStatus& contact_status,
                         const SplitSolution& s, ContactDynamicsData& data);

///
/// @brief Computes the residual and derivatives of the contact dynamics  
/// constraint and derivatives of it. 
/// @param[in] robot Robot model. 
/// @param[in] contact_status Contact status of this time stage. 
/// @param[in] s Split solution of this time stage.
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in, out] kkt_residual Split KKT residual of this time stage.
///
void linearizeContactDynamics(Robot& robot, const ContactStatus& contact_status, 
                              const SplitSolution& s, 
                              ContactDynamicsData& data,
                              SplitKKTResidual& kkt_residual);

///
/// @brief Condenses the acceleration, contact forces, and Lagrange
/// multipliers. 
/// @param[in] robot Robot model. 
/// @param[in] contact_status Contact status of this time stage. 
/// @param[in] dt Time step of this time stage. 
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
/// @param[in, out] kkt_residual Split KKT residual of this time stage.
///
void condenseContactDynamics(Robot& robot, const ContactStatus& contact_status, 
                             const double dt, ContactDynamicsData& data,
                             SplitKKTMatrix& kkt_matrix, 
                             SplitKKTResidual& kkt_residual);

///
/// @brief Condenses the switching constraint. 
/// @param[in] data Data structure for the contact dynamics.
/// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
/// @param[in, out] sc_residual Residual of the switching constraint. 
///
void condenseContactDynamics(const ContactDynamicsData& data, 
                             SwitchingConstraintJacobian& sc_jacobian,
                             SwitchingConstraintResidual& sc_residual);

///
/// @brief Expands the primal variables, i.e., computes the Newton direction 
/// of the condensed primal variables (acceleration a and the contact forces 
/// f) of this stage.
/// @param[in] data Data structure for the contact dynamics.
/// @param[in, out] d Split direction of this time stage.
/// 
void expandContactDynamicsPrimal(const ContactDynamicsData& data, SplitDirection& d);

///
/// @brief Expands the dual variables, i.e., computes the Newton direction 
/// of the condensed dual variables (Lagrange multipliers) of this stage.
/// @param[in] dt Time step of this time stage. 
/// @param[in] dts Direction of the switching time regarding of this time stage. 
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in] d_next Split direction of the next stage.
/// @param[in, out] d Split direction of this time stage.
/// 
void expandContactDynamicsDual(const double dt, const double dts, 
                               ContactDynamicsData& data,
                               const SplitDirection& d_next, SplitDirection& d);

///
/// @brief Expands the dual variables, i.e., computes the Newton direction 
/// of the condensed dual variables (Lagrange multipliers) of this stage.
/// @param[in] dt Time step of this time stage. 
/// @param[in] dts Direction of the switching time regarding of this time stage. 
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in] sc_jacobian Jacobian of the switching constraint. 
/// @param[in] d_next Split direction of the next stage.
/// @param[in, out] d Split direction of this time stage.
/// 
void expandContactDynamicsDual(const double dt, const double dts, 
                               ContactDynamicsData& data, 
                               const SwitchingConstraintJacobian& sc_jacobian,
                               const SplitDirection& d_next, 
                               SplitDirection& d);

} // namespace robotoc 

#endif // ROBOTOC_CONTACT_DYNAMICS_HPP_ 