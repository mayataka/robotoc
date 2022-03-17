import robotoc
import numpy as np
from a1_simulator import A1Simulator


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id = robot.contact_frames()

step_length = np.array([0.15, 0, 0]) 
# step_length = np.array([-0.1, 0, 0]) 
# step_length = np.array([0, 0.1, 0]) 
# step_length = np.array([0.1, -0.1, 0]) 

step_height = 0.1
swing_time = 0.25
stance_time = 0
stance_time = 0.05
swing_start_time = 0.5

vcom_cmd = step_length / swing_time
yaw_cmd = 0
# yaw_cmd = np.pi / 12

T = 0.5
N = 18
max_steps = 3
nthreads = 4
mpc = robotoc.MPCTrotting(robot, T, N, max_steps, nthreads)

planner = robotoc.TrottingFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (yaw_cmd*swing_time), (stance_time > 0.))
mpc.set_gait_pattern(planner, step_height, swing_time, stance_time, swing_start_time)

q = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3, 
              0.0,  0.67, -1.3])
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
sim_end_time = 50.0
sim = A1Simulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=False)
# sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=False)
# sim.run_simulation(mpc, q, v, verbose=False, record=True, record_name='a1_trotting.mp4')