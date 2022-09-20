import robotoc
from robotoc_sim import MPCSimulation, CameraSettings, TerrainSettings
from a1_simulator import A1Simulator
import numpy as np


model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../a1_description/urdf/a1.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.point_contacts = [robotoc.ContactModelInfo('FL_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('RL_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('FR_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('RR_foot', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

step_length = np.array([0.25, 0, 0]) 
step_yaw = 0.0

step_height = 0.1
swing_time = 0.25
stance_time = 0.05
swing_start_time = 0.5

vcom_cmd = 0.25 * step_length / (swing_time+stance_time)
yaw_rate_cmd = step_yaw / (swing_time+stance_time)

T = 0.5
N = 18
nthreads = 4
mpc = robotoc.MPCCrawl(robot, T, N, nthreads)

planner = robotoc.CrawlFootStepPlanner(robot)
planner.set_gait_pattern(step_length, step_yaw, (stance_time > 0.))
# planner.set_raibert_gait_pattern(vcom_cmd, yaw_rate_cmd, swing_time, stance_time, gain=0.7)
mpc.set_gait_pattern(planner, step_height, swing_time, stance_time, swing_start_time)

t0 = 0.0
q0 = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3])
q0[0] -= 2.5
v0 = np.zeros(robot.dimv())
option_init = robotoc.SolverOptions()
option_init.max_iter = 10
mpc.init(t0, q0, v0, option_init)

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 2 # MPC iterations
mpc.set_solver_options(option_mpc)

time_step = 0.0025 # 400 Hz MPC
a1_simulator = A1Simulator(urdf_path=model_info.urdf_path, time_step=time_step)
terrain_settings = TerrainSettings(from_urdf=True)
a1_simulator.set_terrain_settings(terrain_settings)
camera_settings = CameraSettings(camera_distance=3.0, camera_yaw=15, camera_pitch=8.0, 
                                 camera_target_pos=q0[0:3]+np.array([0.0, 1.5, 0.2]))
a1_simulator.set_camera_settings(camera_settings=camera_settings)

simulation_time = 10.0
log = False
record = False
q0[2] += 0.04 # to avoid penetration at the initial configuraion
simulation = MPCSimulation(simulator=a1_simulator)
simulation.run(mpc=mpc, t0=t0, q0=q0, simulation_time=simulation_time, 
               feedback_delay=True, verbose=False, 
               record=record, log=log, name='a1_crawl_terrain')

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
