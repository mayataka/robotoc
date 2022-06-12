import robotoc
import numpy as np
from a1_simulator import A1Simulator


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

step_length = np.array([0.2, 0, 0]) 
yaw_cmd = 0.0

step_height = 0.1
swing_time = 0.25
# stance_time = 0
stance_time = 0.05
swing_start_time = 0.5

vcom_cmd = step_length / swing_time
yaw_rate_cmd = yaw_cmd / swing_time

T = 0.5
N = 18
max_steps = 3
nthreads = 4
mpc = robotoc.MPCCrawl(robot, T, N, max_steps, nthreads)

planner = robotoc.CrawlFootStepPlanner(robot)
raibert_gain = 0.2
planner.set_gait_pattern(vcom_cmd, yaw_rate_cmd, swing_time, swing_time+stance_time, raibert_gain)
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
option_mpc.max_iter = 2 # MPC iterations
mpc.set_solver_options(option_mpc)

sim_time_step = 0.0025 # 400 Hz MPC
sim_start_time = 0.0
sim_end_time = 10.0
sim = A1Simulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

q0 = q.copy()
q0[2] += 0.1 # to avoid penetration at the initial configuraion
sim.set_camera(3.0, 15, 8, q[0:3]+np.array([0.0, 1.5, 0.2]))
sim.run_simulation(mpc, q0, v, feedback_delay=True, verbose=False, terrain=True, record=False)
# sim.run_simulation(mpc, q, v, verbose=False, record=True, terrain=True, record_name='a1_crawl_terrain.mp4')