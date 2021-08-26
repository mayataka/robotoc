import idocp
import numpy as np
from anymal_simulator import ANYmalSimulator


LF_foot_id = 12
LH_foot_id = 22
RF_foot_id = 32
RH_foot_id = 42
contact_frames = [LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id] 
path_to_urdf = '../anymal_b_simple_description/urdf/anymal.urdf'
baumgarte_time_step = 0.05
robot = idocp.Robot(path_to_urdf, idocp.BaseJointType.FloatingBase, 
                    contact_frames, baumgarte_time_step)

dt = 0.02
step_length = 0.25
step_height = 0.15
period_swing = 0.5
t0 = 0.5

cost = idocp.CostFunction()
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
config_cost = idocp.ConfigurationSpaceCost(robot)
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
q0_3d_LF = robot.frame_position(LF_foot_id)
q0_3d_LH = robot.frame_position(LH_foot_id)
q0_3d_RF = robot.frame_position(RF_foot_id)
q0_3d_RH = robot.frame_position(RH_foot_id)
LF_t0 = t0 + 3 * period_swing 
LH_t0 = t0 + 2 * period_swing 
RF_t0 = t0 + period_swing
RH_t0 = t0
LF_foot_ref = idocp.PeriodicFootTrackRef2(q0_3d_LF, step_length, step_height, 
                                          LF_t0, period_swing, 3*period_swing, False)
LH_foot_ref = idocp.PeriodicFootTrackRef2(q0_3d_LH, step_length, step_height, 
                                          LH_t0, period_swing, 3*period_swing, False)
RF_foot_ref = idocp.PeriodicFootTrackRef2(q0_3d_RF, step_length, step_height, 
                                          RF_t0, period_swing, 3*period_swing, True)
RH_foot_ref = idocp.PeriodicFootTrackRef2(q0_3d_RH, step_length, step_height, 
                                          RH_t0, period_swing, 3*period_swing, True)
LF_cost = idocp.TimeVaryingTaskSpace3DCost(robot, LF_foot_id, LF_foot_ref)
LH_cost = idocp.TimeVaryingTaskSpace3DCost(robot, LH_foot_id, LH_foot_ref)
RF_cost = idocp.TimeVaryingTaskSpace3DCost(robot, RF_foot_id, RF_foot_ref)
RH_cost = idocp.TimeVaryingTaskSpace3DCost(robot, RH_foot_id, RH_foot_ref)
foot_track_weight = np.full(3, 1.0e03)
LF_cost.set_q_weight(foot_track_weight)
LH_cost.set_q_weight(foot_track_weight)
RF_cost.set_q_weight(foot_track_weight)
RH_cost.set_q_weight(foot_track_weight)
cost.push_back(LF_cost)
cost.push_back(LH_cost)
cost.push_back(RF_cost)
cost.push_back(RH_cost)

com_ref0 = (q0_3d_LF + q0_3d_LH + q0_3d_RF + q0_3d_RH) / 4
com_ref0[2] = robot.com()[2]
v_com_ref = np.zeros(3)
v_com_ref[0] = 0.25 * step_length / period_swing
com_ref = idocp.PeriodicCoMRef2(com_ref0, v_com_ref, t0, 2*period_swing, 0., True)
com_cost = idocp.TimeVaryingCoMCost(robot, com_ref)
com_cost.set_q_weight(np.full(3, 1.0e04))
cost.push_back(com_cost)

constraints           = idocp.Constraints()
joint_position_lower  = idocp.JointPositionLowerLimit(robot)
joint_position_upper  = idocp.JointPositionUpperLimit(robot)
joint_velocity_lower  = idocp.JointVelocityLowerLimit(robot)
joint_velocity_upper  = idocp.JointVelocityUpperLimit(robot)
joint_torques_lower   = idocp.JointTorquesLowerLimit(robot)
joint_torques_upper   = idocp.JointTorquesUpperLimit(robot)
mu = 0.7
friction_cone         = idocp.FrictionCone(robot, mu)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)
constraints.set_barrier(1.0e-01)

T = 0.5
N = 20
max_steps = 3

nthreads = 4
mpc = idocp.MPCQuadrupedalWalking(robot, cost, constraints, T, N, 
                                  max_steps, nthreads)
mpc.set_gait_pattern(step_length, step_height, period_swing, t0)
q = q_standing
v = np.zeros(robot.dimv())
t = 0.0
mpc.init(t, q, v, 5)

sim_time_step = 0.00333
sim_start_time = 0.0
sim_end_time = 10.0
sim = ANYmalSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.5, 0., 0.]))
sim.run_simulation(mpc, q, v, num_mpc_iteration=5, verbose=False, record=False)
# sim.run_simulation(mpc, q, v, num_mpc_iteration=5, verbose=False, record=True, record_name='anymal_walking.mp4')