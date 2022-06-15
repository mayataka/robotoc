import robotoc
import numpy as np
from icub_simulator import iCubSimulator


path_to_urdf = '../icub_description/urdf/icub_lower_half.urdf'
contact_frames = ['l_sole', 'r_sole']
contact_types = [robotoc.ContactType.SurfaceContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

knee_angle = np.pi / 6
step_length = np.array([0.225, 0, 0]) 
yaw_cmd = np.pi / 60

step_height = 0.1
swing_time = 0.7
double_support_time = 0.0
# double_support_time = 0.05
swing_start_time = 0.5

vcom_cmd = step_length / swing_time
yaw_rate_cmd = yaw_cmd / swing_time

T = 0.7
N = 20
max_steps = 3
nthreads = 4
mpc = robotoc.MPCBipedWalk(robot, T, N, max_steps, nthreads)

planner = robotoc.BipedWalkFootStepPlanner(robot)
planner.set_gait_pattern(step_length, yaw_cmd, (double_support_time > 0.))
# raibert_gain = 0.5
# planner.set_gait_pattern(vcom_cmd, yaw_rate_cmd, swing_time, swing_time+double_support_time, raibert_gain)
mpc.set_gait_pattern(planner, step_height, swing_time, double_support_time, swing_start_time)

# wrench cone
mu = 0.4
X = 0.05
Y = 0.025
mpc.get_wrench_cone_handle().set_friction_coefficient(mu=mu)
mpc.get_wrench_cone_handle().set_rectangular(X=X, Y=Y)
mpc.get_impulse_wrench_cone_handle().set_friction_coefficient(mu=mu)
mpc.get_impulse_wrench_cone_handle().set_rectangular(X=X, Y=Y)

q = np.array([0, 0, 0, 0, 0, 0, 1,
              0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0,  # left leg
              0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0]) # right leg
robot.forward_kinematics(q)
q[2] = - 0.5 * (robot.frame_position('l_sole')[2] + robot.frame_position('r_sole')[2]) 
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
sim_end_time = 20.0
sim = iCubSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

log = False
record = False

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.5, 1.2, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, 
                   record=record, log=log, sim_name='icub_walk')

if record:
    robotoc.utils.adjust_video_duration('icub_walk.mp4', 
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
    plot_mpc.plot(t_log, v_com_log, w_com_log, fig_name='icub_walk_com_vel')
