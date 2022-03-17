import robotoc
import numpy as np
from icub_simulator import iCubSimulator


path_to_urdf = '../icub_description/urdf/icub_lower_half.urdf'
contact_frames = ['l_sole', 'r_sole']
contact_types = [robotoc.ContactType.SurfaceContact for i in contact_frames]
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
L_foot_id, R_foot_id = robot.contact_frames()

knee_angle = np.pi / 6
step_length = np.array([0.2, 0, 0]) 
step_height = 0.1
swing_time = 0.5
double_support_time = 0.0
# double_support_time = 0.05
swing_start_time = 0.5

vcom_cmd = step_length / swing_time
yaw_cmd = 0

T = 0.5
N = 20
max_steps = 3
nthreads = 4
mpc = robotoc.MPCWalking(robot, T, N, max_steps, nthreads)

planner = robotoc.WalkingFootStepPlanner(robot)
planner.set_gait_pattern(step_length, (yaw_cmd*swing_time), (double_support_time > 0.))
mpc.set_gait_pattern(planner, step_height, swing_time, double_support_time, swing_start_time)

mpc.get_wrench_cone_handle().set_friction_coefficient(mu=0.3)
mpc.get_wrench_cone_handle().set_rectangular(X=0.05, Y=0.025)
mpc.get_impulse_wrench_cone_handle().set_friction_coefficient(mu=0.3)
mpc.get_impulse_wrench_cone_handle().set_rectangular(X=0.05, Y=0.025)

q = np.array([0, 0, 0, 0, 0, 0, 1,
              0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0,  # left leg
              0.5*knee_angle, 0, 0, -knee_angle, 0.5*knee_angle, 0]) # right leg
robot.forward_kinematics(q)
q[2] = - 0.5 * (robot.frame_position(L_foot_id)[2] + robot.frame_position(R_foot_id)[2]) 
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
sim_end_time = 5.0
sim = iCubSimulator(path_to_urdf, sim_time_step, sim_start_time, sim_end_time)

sim.set_camera(2.0, 45, -10, q[0:3]+np.array([0.1, 0.5, 0.]))
sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=True, record=False)
# sim.run_simulation(mpc, q, v, feedback_delay=True, verbose=False, record=False)
