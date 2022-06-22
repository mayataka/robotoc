import robotoc
import numpy as np
from anymal_simulator import ANYmalSimulator


path_to_urdf = '../anymal_b_simple_description/urdf/anymal.urdf'
contact_frames = ['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'] 
contact_types = [robotoc.ContactType.PointContact for i in range(4)]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

step_length = np.array([0.15, 0, 0]) 
step_yaw = 0

swing_height = 0.1
swing_time = 0.25
# stance_time = 0.05
stance_time = 0.0
swing_start_time = 0.5

vcom_cmd = 0.5 * step_length / (swing_time+stance_time)
yaw_rate_cmd = step_yaw / (swing_time+stance_time)

T = 0.5
N = 18
max_steps = 3
nthreads = 4
mpc = robotoc.MPCTrot(robot, T, N, max_steps, nthreads)

planner = robotoc.TrotFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (step_yaw*swing_time), (stance_time > 0.))
mpc.set_gait_pattern(planner, swing_height, swing_time, stance_time, swing_start_time)

q = np.array([0, 0, 0.4842, 0, 0, 0, 1, 
              -0.1,  0.7, -1.0, 
              -0.1, -0.7,  1.0, 
               0.1,  0.7, -1.0, 
               0.1, -0.7,  1.0])
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
sim_end_time = 10.0
sim = ANYmalSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

log = False
record = False

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, 
                   record=record, log=log, sim_name='anymal_trot')

if record:
    sim.disconnect()
    robotoc.utils.adjust_video_duration(sim.sim_name+'.mp4', 
                                        desired_duration_sec=(sim_end_time-sim_start_time))

if log:
    q_log = np.genfromtxt(sim.q_log)
    v_log = np.genfromtxt(sim.v_log)
    t_log = np.genfromtxt(sim.t_log)
    sim_steps = t_log.shape[0]

    vcom_log = []
    wcom_log = []
    vcom_cmd_log = []
    yaw_rate_cmd_log = []
    for i in range(sim_steps):
        R = robotoc.utils.rotation_matrix(q_log[i][3:7])
        robot.forward_kinematics(q_log[i], v_log[i])
        vcom_log.append(R.T@robot.com_velocity()) # robot.com_velocity() is expressed in the world coordinate
        wcom_log.append(v_log[i][3:6])
        vcom_cmd_log.append(vcom_cmd)
        yaw_rate_cmd_log.append(yaw_rate_cmd)

    plot_mpc = robotoc.utils.PlotCoMVelocity()
    plot_mpc.plot(t_log, vcom_log, wcom_log, vcom_cmd_log, yaw_rate_cmd_log, 
                  fig_name=sim.sim_name+'_com_vel')
