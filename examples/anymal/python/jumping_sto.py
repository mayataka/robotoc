import robotoc
import numpy as np
import math


LF_foot_id = 12
LH_foot_id = 22
RF_foot_id = 32
RH_foot_id = 42
contact_frames = [LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id] 
path_to_urdf = '../anymal_b_simple_description/urdf/anymal.urdf'
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, baumgarte_time_step)

dt = 0.02
jump_length = 0.8
flying_up_time = 0.15
flying_down_time = flying_up_time
flying_time = flying_up_time + flying_down_time
ground_time = 0.7
t0 = 0.

# Create the cost function
cost = robotoc.CostFunction()
q_standing = np.array([0., 0., 0.4792, 0., 0., 0., 1.0, 
                       -0.1,  0.7, -1.0, 
                       -0.1, -0.7,  1.0, 
                        0.1,  0.7, -1.0, 
                        0.1, -0.7,  1.0])
q_ref = q_standing.copy()
q_ref[0] += jump_length
q_weight = np.array([1.0, 0., 0., 1.0, 1.0, 1.0, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001])
v_weight = np.full(robot.dimv(), 1.0)
u_weight = np.full(robot.dimu(), 1.0e-06)
a_weight = np.full(robot.dimv(), 1.0e-06)
qi_weight = np.array([0., 0., 0., 100., 100., 100., 
                      0.1, 0.1, 0.1, 
                      0.1, 0.1, 0.1,
                      0.1, 0.1, 0.1,
                      0.1, 0.1, 0.1])
vi_weight = np.full(robot.dimv(), 1.0)
dvi_weight = np.full(robot.dimv(), 1.0e-06)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_ref)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_dvi_weight(dvi_weight)
config_cost.set_a_weight(a_weight)
cost.push_back(config_cost)

robot.forward_kinematics(q_standing)
q0_3d_LF = robot.frame_position(LF_foot_id)
q0_3d_LH = robot.frame_position(LH_foot_id)
q0_3d_RF = robot.frame_position(RF_foot_id)
q0_3d_RH = robot.frame_position(RH_foot_id)

# Create the constraints
constraints           = robotoc.Constraints()
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
mu = 0.7
friction_cone         = robotoc.FrictionCone(robot, mu)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)
constraints.set_barrier(1.0e-03)

# Create the contact sequence
max_num_impulses = 1
contact_sequence = robotoc.ContactSequence(robot, max_num_impulses)

contact_points = [q0_3d_LF, q0_3d_LH, q0_3d_RF, q0_3d_RH]
contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts([0, 1, 2, 3])
contact_status_standing.set_contact_points(contact_points)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_flying = robot.create_contact_status()
contact_sequence.push_back(contact_status_flying, t0+ground_time-0.3, sto=True)

contact_points[0][0] += jump_length
contact_points[1][0] += jump_length
contact_points[2][0] += jump_length
contact_points[3][0] += jump_length
contact_status_standing.set_contact_points(contact_points)
contact_sequence.push_back(contact_status_standing, t0+ground_time+flying_time-0.1, sto=True)

# you can check the contact sequence via 
# print(contact_sequence)

T = t0 + flying_time + 2*ground_time
N = math.floor(T/dt) 
ocp_solver = robotoc.OCPSolver(robot, contact_sequence, cost, constraints, 
                               T, N, nthreads=4)

t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

ocp_solver.set_discretization_method(robotoc.DiscretizationMethod.PhaseBased) 
ocp_solver.mesh_refinement(t)
ocp_solver.set_STO_constraints(1.0e-03, 0.95, 0.1, 0.1, 0.65)

ocp_solver.init_constraints(t)

reg_type = robotoc.STORegularizationType.Quad
reg_w = 0.1
ocp_solver.set_STO_regularization(reg_type, reg_w)

num_iteration = 50
robotoc.utils.benchmark.convergence_sto(ocp_solver, t, q, v, num_iteration)
# num_iteration = 1000
# robotoc.utils.benchmark.cpu_time(ocp_solver, t, q, v, num_iteration)

viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='gepetto')
viewer.set_contact_info(contact_frames, mu)
viewer.display(dt, ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))