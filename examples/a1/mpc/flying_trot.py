import robotoc
import numpy as np
from a1_simulator import A1Simulator


# cmd_type = 'forward'
# cmd_type = 'backward'
# cmd_type = 'side'
# cmd_type = 'curve'
cmd_type = 'rotation'

if cmd_type == 'forward':
    step_length = np.array([0.30, 0.0, 0.0]) 
    yaw_cmd = 0.
elif cmd_type == 'backward':
    step_length = np.array([-0.1, 0.0, 0.0]) 
    yaw_cmd = 0.
elif cmd_type == 'side':
    step_length = np.array([0.0, 0.15, 0.0]) 
    yaw_cmd = 0.
elif cmd_type == 'curve':
    step_length = np.array([0.1, 0.0, 0.0]) 
    yaw_cmd = np.pi / 12.
elif cmd_type == 'rotation':
    step_length = np.array([0.0, 0.0, 0.0]) 
    yaw_cmd = np.pi / 6.


path_to_urdf = '../a1_description/urdf/a1.urdf'
contact_frames = ['FL_foot', 'RL_foot', 'FR_foot', 'RR_foot'] 
contact_types = [robotoc.ContactType.PointContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

step_height = 0.1
stance_time = 0.15
flying_time = 0.1
swing_start_time = 0.5

T = 0.5
N = 18
max_steps = 3
nthreads = 4
mpc = robotoc.MPCFlyingTrot(robot, T, N, max_steps, nthreads)


swing_time = stance_time + flying_time
vcom_cmd = step_length / swing_time
yaw_rate_cmd = yaw_cmd / swing_time

planner = robotoc.FlyingTrotFootStepPlanner(robot)
# planner.set_gait_pattern(step_length, yaw_cmd)
raibert_gain = 0.2
planner.set_gait_pattern(vcom_cmd, yaw_rate_cmd, swing_time, stance_time, raibert_gain)
mpc.set_gait_pattern(planner, step_height, flying_time, stance_time, swing_start_time)

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
# sim_end_time = 5.0

sim = A1Simulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
# sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=True)
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=False)
# sim.run_simulation(mpc, q, v, verbose=False, record=True, record_name='rotation.mp4')