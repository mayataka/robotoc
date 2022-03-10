import robotoc
import numpy as np
from icub_simulator import iCubSimulator


L_foot_id = 26
R_foot_id = 46
contact_frames = [L_foot_id, R_foot_id]
contact_types = [robotoc.ContactType.SurfaceContact for i in contact_frames]
path_to_urdf = '../icub_description/urdf/icub_lower_half.urdf'
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

dt = 0.02
knee_angle = np.pi / 6
step_length = np.array([0.2, 0, 0]) 
step_height = 0.1
swing_time = 0.5
double_support_time = 0.0
# double_support_time = 0.05
initial_lift_time = 0.5

vcom_cmd = step_length / swing_time
yaw_cmd = 0

cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0, 0, 0, 0, 1,
                       0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0,  # left leg
                       0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0]) # right leg
robot.forward_kinematics(q_standing)
height = - 0.5 * (robot.frame_position(L_foot_id)[2] + robot.frame_position(R_foot_id)[2]) 
q_standing[2] = height

cost = robotoc.CostFunction()
q_ref = q_standing.copy()
q_weight = np.array([0, 0, 0, 1000, 1000, 1000, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001])
qf_weight = q_weight
v_weight = np.full(robot.dimv(), 1.0)
u_weight = np.full(robot.dimu(), 1.0e-02)
qi_weight = np.array([0, 0, 0, 1000, 1000, 1000, 
                      1, 1, 1, 1, 1, 1, 
                      1, 1, 1, 1, 1, 1])
vi_weight = np.full(robot.dimv(), 1)
dvi_weight = np.full(robot.dimv(), 1e-02)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_u_weight(u_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_dvi_weight(dvi_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
x3d0_L = robot.frame_position(L_foot_id)
x3d0_R = robot.frame_position(R_foot_id)

L_t0 = initial_lift_time + swing_time + double_support_time
R_t0 = initial_lift_time
L_foot_ref = robotoc.PeriodicFootTrackRef(x3d0_L, step_length, step_height, 
                                          L_t0, swing_time, 
                                          swing_time+double_support_time, False)
R_foot_ref = robotoc.PeriodicFootTrackRef(x3d0_R, step_length, step_height, 
                                          R_t0, swing_time, 
                                          swing_time+double_support_time, True)
L_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, L_foot_id, L_foot_ref)
R_cost = robotoc.TimeVaryingTaskSpace3DCost(robot, R_foot_id, R_foot_ref)
foot_track_weight = np.full(3, 1.0e03)
L_cost.set_x3d_weight(foot_track_weight)
R_cost.set_x3d_weight(foot_track_weight)
cost.push_back(L_cost)
cost.push_back(R_cost)


com_ref0 = robot.com()
vcom_ref = 0.5 * step_length / swing_time
com_ref = robotoc.PeriodicCoMRef(com_ref0, vcom_ref, 
                                 initial_lift_time, swing_time, double_support_time, True)
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
mu = 0.4
wrench_friction_cone  = robotoc.WrenchFrictionCone(robot, mu, 0.05, 0.025)
impulse_wrench_friction_cone  = robotoc.ImpulseWrenchFrictionCone(robot, mu, 0.05, 0.025)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(wrench_friction_cone)
constraints.push_back(impulse_wrench_friction_cone)

T = 0.5
N = 18
max_steps = 3
ocp = robotoc.OCP(robot, cost, constraints, T, N, max_steps)

nthreads = 4
mpc = robotoc.MPCWalking(ocp, nthreads)
mpc.set_gait_pattern(vcom_cmd, yaw_cmd, swing_time, double_support_time, initial_lift_time)

q = q_standing
v = np.zeros(robot.dimv())
t = 0.0
option_init = robotoc.SolverOptions()
option_init.max_iter = 200
mpc.init(t, q, v, option_init)

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 1 # MPC iterations
mpc.set_solver_options(option_mpc)

sim_time_step = 0.0025 # 400 Hz MPC
sim_start_time = 0.0
sim_end_time = 5.0
sim = iCubSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=False)