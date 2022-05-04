import robotoc
import numpy as np
from a1_simulator import A1Simulator


jump_type = 'longitudinal'
# Jump_type = 'lateral'
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
    jump_length = [0.1, 0.0, 0]
    jump_yaw = np.pi / 6


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id = robot.contact_frames()


T = 0.8
N = 18
max_steps = 1
nthreads = 4
mpc = robotoc.MPCJumping(robot, T, N, max_steps, nthreads)

planner = robotoc.JumpingFootStepPlanner(robot)
planner.set_jump_pattern(jump_length, jump_yaw)
mpc.set_jump_pattern(planner, flying_time=0.3, min_flying_time=0.2, 
                     ground_time=0.3, min_ground_time=0.2)

q = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3])
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

sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=False)
# sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=True, record_name=jump_type+'.mp4')