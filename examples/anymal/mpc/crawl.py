import robotoc
from robotoc_sim import MPCSimulation, CameraSettings
from anymal_simulator import ANYmalSimulator
import numpy as np


model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../anymal_b_simple_description/urdf/anymal.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.point_contacts = [robotoc.ContactModelInfo('LF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('LH_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RH_FOOT', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

step_length = np.array([0.15, 0, 0]) 
step_yaw = 0

swing_height = 0.1
swing_time = 0.25
stance_time = 0.05
stance_time = 0.0
swing_start_time = 0.5

vcom_cmd = 0.25 * step_length / (swing_time+stance_time)
yaw_rate_cmd = step_yaw / (swing_time+stance_time)

T = 0.5
N = 20
mpc = robotoc.MPCCrawl(robot, T, N)

planner = robotoc.CrawlFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (step_yaw*swing_time), (stance_time > 0.))
mpc.set_gait_pattern(planner, swing_height, swing_time, stance_time, swing_start_time)

t0 = 0.0
q0 = np.array([0, 0, 0.4842, 0, 0, 0, 1, 
               -0.1,  0.7, -1.0, 
               -0.1, -0.7,  1.0, 
                0.1,  0.7, -1.0, 
                0.1, -0.7,  1.0])
v0 = np.zeros(robot.dimv())

option_init = robotoc.SolverOptions()
option_init.max_iter = 10
option_init.nthreads = 4
mpc.init(t0, q0, v0, option_init)

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 1 # MPC iterations
option_mpc.nthreads = 4
mpc.set_solver_options(option_mpc)

time_step = 0.0025 # 400 Hz MPC
anymal_simulator = ANYmalSimulator(urdf_path=model_info.urdf_path, time_step=time_step)
camera_settings = CameraSettings(camera_distance=2.0, camera_yaw=45, camera_pitch=-10.0, 
                                 camera_target_pos=q0[0:3]+np.array([0.1, 0.5, 0.]))
anymal_simulator.set_camera_settings(camera_settings=camera_settings)

simulation_time = 10.0
log = False
record = False
simulation = MPCSimulation(simulator=anymal_simulator)
simulation.run(mpc=mpc, t0=t0, q0=q0, simulation_time=simulation_time, 
               feedback_delay=True, verbose=False, 
               record=record, log=log, name='anymal_crawl')

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
