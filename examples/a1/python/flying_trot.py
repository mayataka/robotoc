import robotoc
import numpy as np
import math


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in range(4)]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

dt = 0.02
step_length = np.array([0.15, 0, 0])
step_height = 0.1
stance_time = 0.15
flying_time = 0.1
t0 = stance_time
cycle = 5

# Create the cost function
cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3])
q_weight = np.array([0, 0, 0, 10000, 10000, 10000, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001])
v_weight = np.array([100, 100, 100, 100, 100, 100, 
                     1, 1, 1, 
                     1, 1, 1,
                     1, 1, 1,
                     1, 1, 1])
u_weight = np.full(robot.dimu(), 1.0e-01)
q_weight_impulse = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                      100, 100, 100, 
                      100, 100, 100,
                      100, 100, 100,
                      100, 100, 100])
v_weight_impulse = np.full(robot.dimv(), 100)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_q_weight_terminal(q_weight)
config_cost.set_q_weight_impulse(q_weight_impulse)
config_cost.set_v_weight(v_weight)
config_cost.set_v_weight_terminal(v_weight)
config_cost.set_v_weight_impulse(v_weight_impulse)
config_cost.set_u_weight(u_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
x3d0_LF = robot.frame_position('FL_foot')
x3d0_LH = robot.frame_position('RL_foot')
x3d0_RF = robot.frame_position('FR_foot')
x3d0_RH = robot.frame_position('RR_foot')
LF_t0 = t0 + stance_time 
LH_t0 = t0 - flying_time
RF_t0 = t0 - flying_time
RH_t0 = t0 + stance_time
LF_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_LF, step_length, step_height, 
                                           LF_t0, stance_time+2.*flying_time, 
                                           stance_time, False)
LH_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_LH, step_length, step_height, 
                                           LH_t0, stance_time+2.*flying_time, 
                                           stance_time, True)
RF_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_RF, step_length, step_height, 
                                           RF_t0, stance_time+2.*flying_time, 
                                           stance_time, True)
RH_foot_ref = robotoc.PeriodicSwingFootRef(x3d0_RH, step_length, step_height, 
                                           RH_t0, stance_time+2.*flying_time, 
                                           stance_time, False)
LF_cost = robotoc.TaskSpace3DCost(robot, 'FL_foot', LF_foot_ref)
LH_cost = robotoc.TaskSpace3DCost(robot, 'RL_foot', LH_foot_ref)
RF_cost = robotoc.TaskSpace3DCost(robot, 'FR_foot', RF_foot_ref)
RH_cost = robotoc.TaskSpace3DCost(robot, 'RR_foot', RH_foot_ref)
foot_track_weight = np.full(3, 1.0e05)
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
com_cost.set_weight(np.full(3, 1.0e05))
cost.push_back(com_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
mu = 0.5
friction_cone         = robotoc.FrictionCone(robot, mu)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)

# Create the contact sequence
max_num_each_discrete_events = 3*cycle
contact_sequence = robotoc.ContactSequence(robot, max_num_each_discrete_events)

contact_positions = {'FL_foot': x3d0_LF, 'RL_foot': x3d0_LH, 'FR_foot': x3d0_RF, 'RR_foot': x3d0_RH} 
contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts(['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'])
contact_status_standing.set_contact_placements(contact_positions)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_lhrf_swing = robot.create_contact_status()
contact_status_lhrf_swing.activate_contacts(['FL_foot', 'RR_foot'])
contact_status_lhrf_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_lhrf_swing, t0)

contact_status_flying = robot.create_contact_status()
contact_status_flying.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_flying, t0+stance_time)

contact_positions['RL_foot'] += 0.5 * step_length
contact_positions['FR_foot'] += 0.5 * step_length
contact_status_rhlf_swing = robot.create_contact_status()
contact_status_rhlf_swing.activate_contacts(['RL_foot', 'FR_foot'])
contact_status_rhlf_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_rhlf_swing, t0+stance_time+flying_time)

contact_status_flying.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_flying, t0+2*stance_time+flying_time)

contact_positions['FL_foot'] += step_length
contact_positions['RR_foot'] += step_length
contact_status_lhrf_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_lhrf_swing, t0+2*stance_time+2*flying_time)

for j in range(cycle-1):
    i = j + 1
    contact_status_flying.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_flying, t0+(2*i+1)*stance_time+(2*i)*flying_time)

    contact_positions['RL_foot'] += step_length
    contact_positions['FR_foot'] += step_length
    contact_status_rhlf_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_rhlf_swing, t0+(2*i+1)*stance_time+(2*i+1)*flying_time)

    contact_status_flying.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_flying, t0+(2*i+2)*stance_time+(2*i+1)*flying_time)

    contact_positions['FL_foot'] += step_length
    contact_positions['RR_foot'] += step_length
    contact_status_lhrf_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_lhrf_swing, t0+(2*i+2)*stance_time+(2*i+2)*flying_time)

T = t0 + cycle*(2*stance_time+2*flying_time) + stance_time
N = math.floor(T/dt) 
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  contact_sequence=contact_sequence, T=T, N=N)
solver_options = robotoc.SolverOptions()
ocp_solver = robotoc.OCPSolver(ocp=ocp, solver_options=solver_options, nthreads=4)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

ocp_solver.init_constraints(t)
print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
ocp_solver.solve(t, q, v)
print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))
print(ocp_solver.get_solver_statistics())

# num_iteration = 1000
# robotoc.utils.benchmark.cpu_time(ocp_solver, t, q, v, num_iteration)

viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='gepetto')
viewer.set_contact_info(robot.contact_frames(), mu)
discretization = ocp_solver.get_time_discretization()
viewer.display(discretization.time_steps(), ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))