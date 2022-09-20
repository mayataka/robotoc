import robotoc
from robotoc_sim import MPCSimulation, CameraSettings
from a1_simulator import A1Simulator
import numpy as np


jump_type = 'longitudinal'
# jump_type = 'lateral'
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

model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../a1_description/urdf/a1.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.point_contacts = [robotoc.ContactModelInfo('FL_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('RL_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('FR_foot', baumgarte_time_step),
                             robotoc.ContactModelInfo('RR_foot', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

T = 0.8
N = 18
nthreads = 4
mpc = robotoc.MPCJump(robot, T, N, nthreads)

planner = robotoc.JumpFootStepPlanner(robot)
planner.set_jump_pattern(jump_length, jump_yaw)
mpc.set_jump_pattern(planner, flying_time=0.3, min_flying_time=0.2, 
                     ground_time=0.3, min_ground_time=0.2)

t0 = 0.0
q0 = np.array([0, 0, 0.3181, 0, 0, 0, 1, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3, 
               0.0,  0.67, -1.3])
v0 = np.zeros(robot.dimv())
option_init = robotoc.SolverOptions()
option_init.max_iter = 50
option_init.initial_sto_reg_iter = 50
option_init.enable_solution_interpolation = False
mpc.init(t0, q0, v0, option_init, sto=True)  

option_mpc = robotoc.SolverOptions()
option_mpc.max_iter = 2 # MPC iterations
option_mpc.initial_sto_reg_iter = 0
option_mpc.max_dt_mesh = T / N
option_mpc.enable_solution_interpolation = False
mpc.set_solver_options(option_mpc)

time_step = 0.0025 # 400 Hz MPC
a1_simulator = A1Simulator(urdf_path=model_info.urdf_path, time_step=time_step)
if jump_type == 'longitudinal':
    camera_settings = CameraSettings(2.0, 35,  -0, q0[0:3]+np.array([-0.1, 0.5, 0.]))
elif jump_type == 'lateral':
    camera_settings = CameraSettings(2.0, 55,  -0, q0[0:3]+np.array([-0.1, 0.5, 0.]))
elif jump_type == 'back':
    camera_settings = CameraSettings(1.5, 20,  -0, q0[0:3]+np.array([-0.4, 0.3, 0.]))
elif jump_type == 'rotational':
    camera_settings = CameraSettings(2.0, 45, -10, q0[0:3]+np.array([-0.1, 0.5, 0.]))
a1_simulator.set_camera_settings(camera_settings=camera_settings)

simulation_time = 5.0
log = False
record = False
simulation = MPCSimulation(simulator=a1_simulator)
simulation.run(mpc=mpc, t0=t0, q0=q0, simulation_time=simulation_time, 
               feedback_delay=True, verbose=False, 
               record=record, log=log, name='a1_jump')

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
    for i in range(sim_steps):
        R = robotoc.utils.rotation_matrix_from_quaternion(q_log[i][3:7])
        robot.forward_kinematics(q_log[i], v_log[i])
        vcom_log.append(R.T@robot.com_velocity()) # robot.com_velocity() is expressed in the world coordinate
        wcom_log.append(v_log[i][3:6])

    plot_mpc = robotoc.utils.PlotCoMVelocity()
    plot_mpc.plot(t_log, vcom_log, wcom_log, 
                  fig_name=simulation.name+'_com_vel')
