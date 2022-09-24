#ifndef ROBOTOC_IMPULSE_DYNAMICS_HPP_
#define ROBOTOC_IMPULSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"


namespace robotoc {

///
/// @brief Computes the residual in the impulse dynamics constraint. 
/// @param[in] robot Robot model. 
/// @param[in] impulse_status Impulse status of this impulse stage. 
/// @param[in] s Split solution of this impulse stage.
/// @param[in, out] data Data structure for the contact dynamics.
///
void evalImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status,
                         const SplitSolution& s, ContactDynamicsData& data);

///
/// @brief Linearizes the impulse dynamics constraint. 
/// @param[in] robot Robot model. 
/// @param[in] impulse_status Impulse status of this impulse stage. 
/// @param[in] s Split solution of this impulse stage.
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
///
void linearizeImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status, 
                              const SplitSolution& s, ContactDynamicsData& data, 
                              SplitKKTResidual& kkt_residual);

///
/// @brief Condenses the inverse dynamics constraint. 
/// @param[in] robot Robot model. 
/// @param[in] impulse_status Impulse status of this impulse stage. 
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
/// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
///
void condenseImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status,
                             ContactDynamicsData& data, 
                             SplitKKTMatrix& kkt_matrix, 
                             SplitKKTResidual& kkt_residual);

///
/// @brief Expands the primal variables, i.e., computes the Newton direction 
/// of the condensed primal variables (impulse change in the velocity dv and 
/// the impulse forces f) of this impulse stage.
/// @param[in] data Data structure for the contact dynamics.
/// @param[in, out] d Split direction of this impulse stage.
/// 
void expandImpulseDynamicsPrimal(const ContactDynamicsData& data, 
                                 SplitDirection& d);

///
/// @brief Expands the dual variables, i.e., computes the Newton direction 
/// of the condensed dual variables (Lagrange multipliers) of this impulse 
/// stage.
/// @param[in, out] data Data structure for the contact dynamics.
/// @param[in] d_next Split direction of the next stage.
/// @param[in, out] d Split direction of this impulse stage.
/// 
void expandImpulseDynamicsDual(ContactDynamicsData& data, 
                               const SplitDirection& d_next, SplitDirection& d);

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_DYNAMICS_HPP_ 