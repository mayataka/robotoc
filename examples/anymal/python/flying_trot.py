import robotoc
import numpy as np
import math


model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../anymal_b_simple_description/urdf/anymal.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.point_contacts = [robotoc.ContactModelInfo('LF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('LH_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RH_FOOT', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

dt = 0.02
step_length = np.array([0.25, 0, 0])
step_height = 0.1
stance_time = 0.2
flying_time = 0.1
t0 = stance_time
cycle = 3

# Create the cost function
cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.4792, 0, 0, 0, 1, 
                       -0.1,  0.7, -1.0, 
                       -0.1, -0.7,  1.0, 
                        0.1,  0.7, -1.0, 
                        0.1, -0.7,  1.0])
q_weight = np.array([0, 0, 0, 250000, 250000, 250000, 
                     0.0001, 0.0001, 0.0001, 
                     0.0001, 0.0001, 0.0001,
                     0.0001, 0.0001, 0.0001,
                     0.0001, 0.0001, 0.0001])
v_weight = np.array([100, 100, 100, 100, 100, 100, 
                     1, 1, 1, 
                     1, 1, 1,
                     1, 1, 1,
                     1, 1, 1])
u_weight = np.full(robot.dimu(), 1.0e-01)
q_weight_impact = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                      100, 100, 100, 
                      100, 100, 100,
                      100, 100, 100,
                      100, 100, 100])
v_weight_impact = np.full(robot.dimv(), 100)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_q_weight_terminal(q_weight)
config_cost.set_q_weight_impact(q_weight_impact)
config_cost.set_v_weight(v_weight)
config_cost.set_v_weight_terminal(v_weight)
config_cost.set_v_weight_impact(v_weight_impact)
config_cost.set_u_weight(u_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
x3d0_LF = robot.frame_position('LF_FOOT')
x3d0_LH = robot.frame_position('LH_FOOT')
x3d0_RF = robot.frame_position('RF_FOOT')
x3d0_RH = robot.frame_position('RH_FOOT')
LF_t0 = t0 + stance_time 
LH_t0 = t0 - flying_time
RF_t0 = t0 - flying_time
RH_t0 = t0 + stance_time
LF_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_LF, step_length, step_height, 
                                           LF_t0, stance_time+2*flying_time, 
                                           stance_time, False)
LH_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_LH, step_length, step_height, 
                                           LH_t0, stance_time+2*flying_time, 
                                           stance_time, True)
RF_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_RF, step_length, step_height, 
                                           RF_t0, stance_time+2*flying_time, 
                                           stance_time, True)
RH_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_RH, step_length, step_height, 
                                           RH_t0, stance_time+2*flying_time, 
                                           stance_time, False)
LF_cost = robotoc.TaskSpace3DCost(robot, 'LF_FOOT', LF_foot_ref)
LH_cost = robotoc.TaskSpace3DCost(robot, 'LH_FOOT', LH_foot_ref)
RF_cost = robotoc.TaskSpace3DCost(robot, 'RF_FOOT', RF_foot_ref)
RH_cost = robotoc.TaskSpace3DCost(robot, 'RH_FOOT', RH_foot_ref)
foot_track_weight = np.full(3, 1.0e06)
LF_cost.set_weight(foot_track_weight)
LH_cost.set_weight(foot_track_weight)
RF_cost.set_weight(foot_track_weight)
RH_cost.set_weight(foot_track_weight)
cost.push_back(LF_cost)
cost.push_back(LH_cost)
cost.push_back(RF_cost)
cost.push_back(RH_cost)

com_ref0 = robot.com()
vcom_ref = 0.5 * step_length / (stance_time+flying_time)
com_ref = robotoc.PeriodicCoMRef(com_ref0, vcom_ref, t0, stance_time+flying_time, 
                                 0, True)
com_cost = robotoc.CoMCost(robot, com_ref)
com_cost.set_weight(np.full(3, 1.0e06))
cost.push_back(com_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier_param=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
friction_cone         = robotoc.FrictionCone(robot)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)

# Create the contact sequence
contact_sequence = robotoc.ContactSequence(robot)
mu = 0.7
friction_coefficients = {'LF_FOOT': mu, 'LH_FOOT': mu, 'RF_FOOT': mu, 'RH_FOOT': mu} 

contact_positions = {'LF_FOOT': x3d0_LF, 'LH_FOOT': x3d0_LH, 'RF_FOOT': x3d0_RF, 'RH_FOOT': x3d0_RH} 
contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts(['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'])
contact_status_standing.set_contact_placements(contact_positions)
contact_status_standing.set_friction_coefficients(friction_coefficients)
contact_sequence.init(contact_status_standing)

contact_status_lhrf_swing = robot.create_contact_status()
contact_status_lhrf_swing.activate_contacts(['LF_FOOT', 'RH_FOOT'])
contact_status_lhrf_swing.set_contact_placements(contact_positions)
contact_status_lhrf_swing.set_friction_coefficients(friction_coefficients)
contact_sequence.push_back(contact_status_lhrf_swing, t0)

contact_status_flying = robot.create_contact_status()
contact_status_flying.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_flying, t0+stance_time)

contact_positions['LH_FOOT'] += 0.5 * step_length
contact_positions['RF_FOOT'] += 0.5 * step_length
contact_status_rhlf_swing = robot.create_contact_status()
contact_status_rhlf_swing.activate_contacts(['LH_FOOT', 'RF_FOOT'])
contact_status_rhlf_swing.set_contact_placements(contact_positions)
contact_status_rhlf_swing.set_friction_coefficients(friction_coefficients)
contact_sequence.push_back(contact_status_rhlf_swing, t0+stance_time+flying_time)

contact_status_flying.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_flying, t0+2*stance_time+flying_time)

contact_positions['LF_FOOT'] += step_length
contact_positions['RH_FOOT'] += step_length
contact_status_lhrf_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_lhrf_swing, t0+2*stance_time+2*flying_time)

for j in range(cycle-1):
    i = j + 1
    contact_status_flying.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_flying, t0+(2*i+1)*stance_time+(2*i)*flying_time)

    contact_positions['LH_FOOT'] += step_length
    contact_positions['RF_FOOT'] += step_length
    contact_status_rhlf_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_rhlf_swing, t0+(2*i+1)*stance_time+(2*i+1)*flying_time)

    contact_status_flying.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_flying, t0+(2*i+2)*stance_time+(2*i+1)*flying_time)

    contact_positions['LF_FOOT'] += step_length
    contact_positions['RH_FOOT'] += step_length
    contact_status_lhrf_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_lhrf_swing, t0+(2*i+2)*stance_time+(2*i+2)*flying_time)

T = t0 + cycle*(2*stance_time+2*flying_time) + stance_time
N = math.floor(T/dt) 
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  contact_sequence=contact_sequence, T=T, N=N)
solver_options = robotoc.SolverOptions()
solver_options.nthreads = 4
ocp_solver = robotoc.OCPSolver(ocp=ocp, solver_options=solver_options)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.discretize(t)
ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

ocp_solver.init_constraints()
print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
ocp_solver.solve(t, q, v)
print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))
print(ocp_solver.get_solver_statistics())

# num_iteration = 1000
# robotoc.utils.benchmark.cpu_time(ocp_solver, t, q, v, num_iteration)

viewer = robotoc.utils.TrajectoryViewer(model_info=model_info, viewer_type='gepetto')
viewer.set_contact_info(mu=mu)
viewer.display(ocp_solver.get_time_discretization(), 
               ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))