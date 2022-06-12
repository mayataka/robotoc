import robotoc
import numpy as np
import math


path_to_urdf = '../anymal_b_simple_description/urdf/anymal.urdf'
contact_frames = ['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'] 
contact_types = [robotoc.ContactType.PointContact for i in range(4)]
baumgarte_time_step = 0.04
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id = robot.contact_frames()

dt = 0.02
step_length = np.array([0.15, 0, 0])
step_height = 0.1
swing_time = 0.5
double_support_time = 0.04
t0 = double_support_time 
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
qi_weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                      100, 100, 100, 
                      100, 100, 100,
                      100, 100, 100,
                      100, 100, 100])
vi_weight = np.full(robot.dimv(), 100)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_u_weight(u_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
x3d0_LF = robot.frame_position(LF_foot_id)
x3d0_LH = robot.frame_position(LH_foot_id)
x3d0_RF = robot.frame_position(RF_foot_id)
x3d0_RH = robot.frame_position(RH_foot_id)
LF_foot_ref = robotoc.DiscreteTimeSwingFootRef(contact_index=0, swing_height=step_height)
LH_foot_ref = robotoc.DiscreteTimeSwingFootRef(contact_index=1, swing_height=step_height)
RF_foot_ref = robotoc.DiscreteTimeSwingFootRef(contact_index=2, swing_height=step_height)
RH_foot_ref = robotoc.DiscreteTimeSwingFootRef(contact_index=3, swing_height=step_height)
LF_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, LF_foot_id, LF_foot_ref)
LH_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, LH_foot_id, LH_foot_ref)
RF_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, RF_foot_id, RF_foot_ref)
RH_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, RH_foot_id, RH_foot_ref)
foot_track_weight = np.full(3, 1.0e06)
LF_cost.set_x3d_weight(foot_track_weight)
LH_cost.set_x3d_weight(foot_track_weight)
RF_cost.set_x3d_weight(foot_track_weight)
RH_cost.set_x3d_weight(foot_track_weight)
cost.push_back(LF_cost)
cost.push_back(LH_cost)
cost.push_back(RF_cost)
cost.push_back(RH_cost)

com0 = robot.com()
com_pos_to_LF_foot_pos = x3d0_LF - com0
com_pos_to_LH_foot_pos = x3d0_LH - com0
com_pos_to_RF_foot_pos = x3d0_RF - com0
com_pos_to_RH_foot_pos = x3d0_RH - com0
com_pos_to_fee_pos = [com_pos_to_LF_foot_pos, 
                      com_pos_to_LH_foot_pos,
                      com_pos_to_RF_foot_pos, 
                      com_pos_to_RH_foot_pos]
com_ref = robotoc.DiscreteTimeCoMRef(com_pos_to_fee_pos)
com_cost = robotoc.TimeVaryingCoMCost(robot, com_ref)
com_cost.set_com_weight(np.full(3, 1.0e04))
cost.push_back(com_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
mu = 0.7
friction_cone         = robotoc.FrictionCone(robot, mu)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)

# Create the contact sequence
max_num_each_discrete_events = 2*cycle
contact_sequence = robotoc.ContactSequence(robot, max_num_each_discrete_events)

contact_positions = [x3d0_LF, x3d0_LH, x3d0_RF, x3d0_RH]
contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts(['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'])
contact_status_standing.set_contact_placements(contact_positions)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_lhrf_swing = robot.create_contact_status()
contact_status_lhrf_swing.activate_contacts(['LF_FOOT', 'RH_FOOT'])
contact_status_lhrf_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_lhrf_swing, t0, sto=True)

contact_positions[1] += 0.5 * step_length
contact_positions[2] += 0.5 * step_length
contact_status_standing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_standing, t0+swing_time, sto=True)

contact_status_lfrh_swing = robot.create_contact_status()
contact_status_lfrh_swing.activate_contacts(['LH_FOOT', 'RF_FOOT'])
contact_status_lfrh_swing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_lfrh_swing, 
                           t0+swing_time+double_support_time, sto=True)

contact_positions[0] += step_length
contact_positions[3] += step_length
contact_status_standing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_standing, 
                           t0+2*swing_time+double_support_time, sto=True)

for i in range(cycle-1):
    t1 = t0 + (i+1)*(2*swing_time+2*double_support_time)
    contact_status_lhrf_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_lhrf_swing, t1, sto=True)

    contact_positions[1] += step_length
    contact_positions[2] += step_length
    contact_status_standing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_standing, t1+swing_time, sto=True)

    contact_status_lfrh_swing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_lfrh_swing, 
                               t1+swing_time+double_support_time, sto=True)

    contact_positions[0] += step_length
    contact_positions[3] += step_length
    contact_status_standing.set_contact_placements(contact_positions)
    contact_sequence.push_back(contact_status_standing, 
                               t1+2*swing_time+double_support_time, sto=True)

# Sets the swing foot and com refs from the contact sequence
LF_foot_ref.set_swing_foot_ref(contact_sequence)
LH_foot_ref.set_swing_foot_ref(contact_sequence)
RF_foot_ref.set_swing_foot_ref(contact_sequence)
RH_foot_ref.set_swing_foot_ref(contact_sequence)
com_ref.set_com_ref(contact_sequence)

# you can chech the contact sequence as 
# print(contact_sequence)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
min_dt = [0.02] + cycle * [0.2, 0.02, 0.2, 0.02]
sto_constraints = robotoc.STOConstraints(max_num_switches=2*max_num_each_discrete_events, 
                                         min_dt=min_dt,
                                         barrier=1.0e-03, 
                                         fraction_to_boundary_rule=0.995)

T = t0 + cycle*(2*double_support_time+2*swing_time)
N = math.floor(T/dt) 
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  T=T, N=N, max_num_each_discrete_events=max_num_each_discrete_events)
solver_options = robotoc.SolverOptions()
solver_options.max_iter = 200
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
solver_options.initial_sto_reg_iter = 10 
ocp_solver = robotoc.OCPSolver(ocp=ocp, contact_sequence=contact_sequence, 
                               solver_options=solver_options, nthreads=4)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

print('aa')
ocp_solver.mesh_refinement(t)
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