import robotoc
import numpy as np
from anymal_simulator import ANYmalSimulator


contact_frames = ['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'] 
contact_types = [robotoc.ContactType.PointContact for i in range(4)]
path_to_urdf = '../anymal_b_simple_description/urdf/anymal.urdf'
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id = robot.contact_frames()


dt = 0.02
step_length = np.array([0.25, 0, 0]) 
step_height = 0.15
swing_time = 0.25
initial_lift_time = 0.5

vcom_cmd = step_length / swing_time
yaw_cmd = 0

cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.4842, 0, 0, 0, 1, 
                       -0.1,  0.7, -1.0, 
                       -0.1, -0.7,  1.0, 
                        0.1,  0.7, -1.0, 
                        0.1, -0.7,  1.0])
q_weight = np.array([0, 0, 0, 100, 100, 100, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001])
v_weight = np.array([1, 1, 1, 1, 1, 1, 
                     1, 1, 1, 
                     1, 1, 1,
                     1, 1, 1,
                     1, 1, 1])
u_weight = np.full(robot.dimu(), 1.0e-02)
qi_weight = np.array([0, 0, 0, 100, 100, 100, 
                      1, 1, 1, 
                      1, 1, 1,
                      1, 1, 1,
                      1, 1, 1])
vi_weight = np.full(robot.dimv(), 1)
dvi_weight = np.full(robot.dimv(), 1e-03)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_dvi_weight(dvi_weight)
config_cost.set_u_weight(u_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
x3d_LF = robot.frame_position(LF_foot_id)
x3d_LH = robot.frame_position(LH_foot_id)
x3d_RF = robot.frame_position(RF_foot_id)
x3d_RH = robot.frame_position(RH_foot_id)
LF_t0 = initial_lift_time + 3 * swing_time 
LH_t0 = initial_lift_time + 2 * swing_time 
RF_t0 = initial_lift_time + swing_time
RH_t0 = initial_lift_time
LF_foot_ref = robotoc.PeriodicFootTrackRef(x3d_LF, step_length, step_height, 
                                           LF_t0, swing_time, 3*swing_time, False)
LH_foot_ref = robotoc.PeriodicFootTrackRef(x3d_LH, step_length, step_height, 
                                           LH_t0, swing_time, 3*swing_time, False)
RF_foot_ref = robotoc.PeriodicFootTrackRef(x3d_RF, step_length, step_height, 
                                           RF_t0, swing_time, 3*swing_time, True)
RH_foot_ref = robotoc.PeriodicFootTrackRef(x3d_RH, step_length, step_height, 
                                           RH_t0, swing_time, 3*swing_time, True)
LF_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, LF_foot_id, LF_foot_ref)
LH_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, LH_foot_id, LH_foot_ref)
RF_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, RF_foot_id, RF_foot_ref)
RH_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, RH_foot_id, RH_foot_ref)
foot_track_weight = np.full(3, 1.0e03)
LF_cost.set_x3d_weight(foot_track_weight)
LH_cost.set_x3d_weight(foot_track_weight)
RF_cost.set_x3d_weight(foot_track_weight)
RH_cost.set_x3d_weight(foot_track_weight)
cost.push_back(LF_cost)
cost.push_back(LH_cost)
cost.push_back(RF_cost)
cost.push_back(RH_cost)

com_ref0 = robot.com()
vcom_ref = 0.25 * step_length / swing_time
com_ref = robotoc.PeriodicCoMRef(com_ref0, vcom_ref, initial_lift_time, 2*swing_time, 0., True)
com_cost = robotoc.TimeVaryingCoMCost(robot, com_ref)
com_cost.set_com_weight(np.full(3, 1.0e03))
cost.push_back(com_cost)

constraints           = robotoc.Constraints(barrier=1.0e-03)
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

T = 0.5
N = 18
max_steps = 3
ocp = robotoc.OCP(robot, cost, constraints, T, N, max_steps)

planner = robotoc.CrawlingFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (yaw_cmd*swing_time))

nthreads = 4
mpc = robotoc.MPCCrawling(ocp, nthreads)
mpc.set_gait_pattern(planner, swing_time, initial_lift_time)

q = q_standing
v = np.zeros(robot.dimv())
t = 0.0
option_init = robotoc.SolverOptions()
option_init.max_iter = 10
mpc.init(t, q, v, option_init)

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 1 # MPC iterations
mpc.set_solver_options(option_mpc)

sim_time_step = 0.0025 # 400 Hz MPC
sim_start_time = 0.0
sim_end_time = 5.0
sim = ANYmalSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=False)
# sim.run_simulation(mpc, q, v, verbose=False, record=True, record_name='anymal_crawling.mp4')