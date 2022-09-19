import robotoc
from robotoc_sim import MPCSimulation, CameraSettings
from a1_simulator import A1Simulator
import numpy as np


# cmd_type = 'forward'
# cmd_type = 'backward'
# cmd_type = 'side'
# cmd_type = 'curve'
cmd_type = 'rotation'

if cmd_type == 'forward':
    step_length = np.array([0.2, 0.0, 0.0]) 
    step_yaw = 0.
elif cmd_type == 'backward':
    step_length = np.array([-0.1, 0.0, 0.0]) 
    step_yaw = 0.
elif cmd_type == 'side':
    step_length = np.array([0.0, 0.15, 0.0]) 
    step_yaw = 0.
elif cmd_type == 'curve':
    step_length = np.array([0.1, 0.0, 0.0]) 
    step_yaw = np.pi / 30.
elif cmd_type == 'rotation':
    step_length = np.array([0.0, 0.0, 0.0]) 
    step_yaw = np.pi / 20.


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
N = 20
nthreads = 4
mpc = robotoc.MPCFlyingTrot(robot, T, N, nthreads)


vcom_cmd = 0.5 * step_length / (flying_time+stance_time)
yaw_rate_cmd = step_yaw / flying_time

planner = robotoc.FlyingTrotFootStepPlanner(robot)
planner.set_raibert_gait_pattern(vcom_cmd, yaw_rate_cmd, flying_time, stance_time, gain=0.95)
mpc.set_gait_pattern(planner, step_height, flying_time, stance_time, swing_start_time)

t0 = 0.0
q0 = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3])
v0 = np.zeros(robot.dimv())
option_init = robotoc.SolverOptions()
option_init.max_iter = 10
mpc.init(t0, q0, v0, option_init)

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 2 # MPC iterations
mpc.set_solver_options(option_mpc)

time_step = 0.0025 # 400 Hz MPC
a1_simulator = A1Simulator(urdf_path=path_to_urdf, time_step=time_step)
camera_settings = CameraSettings(camera_distance=2.0, camera_yaw=45, camera_pitch=-10.0, 
                                 camera_target_pos=q0[0:3]+np.array([0.1, 0.5, 0.0]))
a1_simulator.set_camera_settings(camera_settings=camera_settings)

simulation_time = 5.0
log = False
record = False
simulation = MPCSimulation(simulator=a1_simulator)
simulation.run(mpc=mpc, t0=t0, q0=q0, simulation_time=simulation_time, 
               feedback_delay=True, verbose=False, 
               record=record, log=log, name='a1_flying_trot')

if record:
    robotoc.utils.adjust_video_duration(simulation.name+'.mp4', 
                                        desired_duration_sec=simulation_time)

if log:
    q_log = np.genfromtxt(simulation.q_log)
    v_log = np.genfromtxt(simulation.v_log)
    t_log = np.genfromtxt(simulation.t_log)
    sim_steps = t_log.shape[0]

    vcom_log = []
    wcom_log = []
    vcom_cmd_log = []
    yaw_rate_cmd_log = []
    for i in range(sim_steps):
        R = robotoc.utils.rotation_matrix_from_quaternion(q_log[i][3:7])
        robot.forward_kinematics(q_log[i], v_log[i])
        vcom_log.append(R.T@robot.com_velocity()) # robot.com_velocity() is expressed in the world coordinate
        wcom_log.append(v_log[i][3:6])
        vcom_cmd_log.append(vcom_cmd)
        yaw_rate_cmd_log.append(yaw_rate_cmd)

    plot_mpc = robotoc.utils.PlotCoMVelocity()
    plot_mpc.plot(t_log, vcom_log, wcom_log, vcom_cmd_log, yaw_rate_cmd_log, 
                  fig_name=simulation.name+'_com_vel')
