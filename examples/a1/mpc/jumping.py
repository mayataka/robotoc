import robotoc
import numpy as np
from a1_simulator import A1Simulator


jump_type = 'longitudinal'
# jump_type = 'lateral'
# jump_type = 'back'
# jump_type = 'rotational'

if jump_type == 'longitudinal':
    jump_length = [0.6, 0, 0]
    jump_yaw = 0
elif jump_type == 'lateral':
    jump_length = [0, 0.4, 0]
    jump_yaw = 0
elif jump_type == 'back':
    jump_length = [-0.3, 0, 0]
    jump_yaw = 0
elif jump_type == 'rotational':
    jump_length = [0.2, 0.0, 0]
    jump_yaw = np.pi / 6


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id = robot.contact_frames()


# Create the cost function
cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3, 
                       0.0,  0.67, -1.3])
robot.normalize_configuration(q_standing)
q_ref = q_standing.copy()
q_ref[0:3] += jump_length

from scipy.spatial.transform import Rotation
rot_mat = [[np.cos(jump_yaw), -np.sin(jump_yaw), 0],
           [np.sin(jump_yaw),  np.cos(jump_yaw), 0],
           [               0,                 0, 1]]
rot_ref = Rotation.from_matrix(rot_mat).as_quat()
q_ref[3:7] = rot_ref
q_weight = np.array([0.0, 100., 100.0, 100., 100., 100., 
                     0.01, 0.01, 0.01, 
                     0.01, 0.01, 0.01,
                     0.01, 0.01, 0.01,
                     0.01, 0.01, 0.01])
v_weight = np.full(robot.dimv(), 1.0)
a_weight = np.full(robot.dimv(), 1.0e-03)
qi_weight = np.array([0., 1000., 1000., 1000., 1000., 1000., 
                      10., 10., 10., 
                      10., 10., 10.,
                      10., 10., 10.,
                      10., 10., 10.])
vi_weight = np.full(robot.dimv(), 10.0)
dvi_weight = np.full(robot.dimv(), 1.0e-03)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_ref)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(qi_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(vi_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_dvi_weight(dvi_weight)
config_cost.set_a_weight(a_weight)
cost.push_back(config_cost)


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


T = 0.8
N = 18
max_steps = 1
sto_cost = robotoc.STOCostFunction()
sto_constraints = robotoc.STOConstraints(max_num_switches=2*max_steps)
ocp = robotoc.OCP(robot, cost, constraints, sto_cost, sto_constraints, 
                  T, N, max_steps)

planner = robotoc.JumpingFootStepPlanner(robot)
planner.set_jump_pattern(jump_length, jump_yaw)

nthreads = 4
mpc = robotoc.MPCJumping(ocp, nthreads)
mpc.set_jump_pattern(planner, flying_time=0.3, min_flying_time=0.2, 
                     ground_time=0.3, min_ground_time=0.2)

q = q_standing
v = np.zeros(robot.dimv())
t = 0.0

option_init = robotoc.SolverOptions()
option_init.max_iter = 50 
option_init.initial_sto_reg_iter = 50
mpc.init(t, q, v, option_init, sto=True)  

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 2 # MPC iterations
option_mpc.initial_sto_reg_iter = 0
option_mpc.max_dt_mesh = T / N
mpc.set_solver_options(option_mpc)

sim_time_step = 0.0025 # 400 Hz MPC
sim_start_time = 0.0
sim_end_time = 1.5
sim = A1Simulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)


if jump_type == 'longitudinal':
    sim.set_camera(2.0, 35, -0, q[0:3]+np.array([-0.1, 0.5, 0.]))
elif jump_type == 'lateral':
    sim.set_camera(2.0, 55, -0, q[0:3]+np.array([-0.1, 0.5, 0.]))
elif jump_type == 'back':
    sim.set_camera(1.5, 20, -0, q[0:3]+np.array([-0.4, 0.3, 0.]))
elif jump_type == 'rotational':
    sim.set_camera(2.0, 45, -10, q[0:3]+np.array([-0.1, 0.5, 0.]))

sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=False)
# sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=True, record_name=jump_type+'.mp4')