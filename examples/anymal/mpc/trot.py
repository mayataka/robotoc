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
swing_height = 0.1
swing_time = 0.25
# stance_time = 0.05
stance_time = 0.0
swing_start_time = 0.5

vcom_cmd = step_length / swing_time
yaw_cmd = 0

T = 0.5
N = 18
max_steps = 3
nthreads = 4
mpc = robotoc.MPCTrot(robot, T, N, max_steps, nthreads)

planner = robotoc.TrotFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (yaw_cmd*swing_time), (stance_time > 0.))
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
    robotoc.utils.adjust_video_duration('anymal_trot.mp4', 
                                        desired_duration_sec=(sim_end_time-sim_start_time))

if log:
    q_log = np.genfromtxt(sim.q_log)
    v_log = np.genfromtxt(sim.v_log)
    t_log = np.genfromtxt(sim.t_log)
    sim_steps = t_log.shape[0]

    from scipy.spatial.transform import Rotation
    v_com_log = []
    w_com_log = []
    for i in range(sim_steps):
        robot.forward_kinematics(q_log[i], v_log[i])
        R = Rotation.from_quat(q_log[i][3:7]).as_matrix()
        v_com_log.append(R.T@robot.com_velocity())
        w_com_log.append(R.T@v_log[i][3:6])

    plot_mpc = robotoc.utils.PlotCoMVelocity()
    plot_mpc.plot(t_log, v_com_log, w_com_log, fig_name='anymal_trot_com_vel')