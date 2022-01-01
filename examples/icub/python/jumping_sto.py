import robotoc
import numpy as np
import math


L_foot_id = 26
R_foot_id = 46
contact_frames = [L_foot_id, R_foot_id]
contact_types = [robotoc.ContactType.SurfaceContact for i in contact_frames]
path_to_urdf = '../icub_description/urdf/icub.urdf'
baumgarte_time_step = 0.05
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)

dt = 0.02
jump_length = 0.5
flying_time = 0.25
ground_time = 0.7
t0 = 0.

cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.592, 0, 0, 1, 0,
                       0.20944, 0.08727, 0, -0.1745, -0.0279, -0.08726, # left leg
                       0.20944, 0.08727, 0, -0.1745, -0.0279, -0.08726, # right leg
                       0, 0, 0, # torso
                       0, 0.35, 0.5, 0.5, 0, 0, 0, # left arm 
                       0, 0.35, 0.5, 0.5, 0, 0, 0]) # right arm 
q_ref = q_standing.copy()
q_ref[0] += jump_length
q_weight = np.array([0, 1, 1, 100, 100, 100, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 1, 1,
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])
qf_weight = q_weight
v_weight = np.full(robot.dimv(), 1.0e-03)
a_weight = np.full(robot.dimv(), 1.0e-05)
qi_weight = 1.0 * q_weight
vi_weight = 1.0 * v_weight
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(qi_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(vi_weight)
config_cost.set_a_weight(a_weight)
cost.push_back(config_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
mu = 0.6
# wrench_friction_cone  = robotoc.WrenchFrictionCone(robot, mu, 0.1, 0.05)
friction_cone  = robotoc.FrictionCone(robot, mu)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
# constraints.push_back(wrench_friction_cone)
constraints.push_back(friction_cone)

# Create the contact sequence
max_num_each_discrete_events = 2
contact_sequence = robotoc.ContactSequence(robot, max_num_each_discrete_events)

robot.forward_kinematics(q_standing)
x3d0_L = robot.frame_placement(L_foot_id)
x3d0_R = robot.frame_placement(R_foot_id)
contact_placements = [x3d0_L, x3d0_R]

contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts([0, 1])
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_flying = robot.create_contact_status()
contact_sequence.push_back(contact_status_flying, t0+ground_time, sto=True)

contact_placements[0].trans = contact_placements[0].trans + np.array([jump_length, 0, 0])
contact_placements[1].trans = contact_placements[1].trans + np.array([jump_length, 0, 0])
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.push_back(contact_status_standing, t0+ground_time+flying_time, sto=True)

contact_sequence.push_back(contact_status_flying, t0+2*ground_time+flying_time, sto=True)

contact_placements[0].trans = contact_placements[0].trans + np.array([jump_length, 0, 0])
contact_placements[1].trans = contact_placements[1].trans + np.array([jump_length, 0, 0])
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.push_back(contact_status_standing, t0+2*ground_time+2*flying_time, sto=True)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
sto_constraints = robotoc.STOConstraints(max_num_switches=2*max_num_each_discrete_events, 
                                         min_dt=[0.6, 0.2, 0.6, 0.2, 0.6],
                                         barrier=1.0e-03, 
                                         fraction_to_boundary_rule=0.995)

T = t0 + 2*flying_time + 3*ground_time
N = math.floor(T/dt) 
# Create the OCP with the STO problem
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  T=T, N=N, max_num_each_discrete_events=max_num_each_discrete_events)
# Create the OCP solver
solver_options = robotoc.SolverOptions()
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
solver_options.max_iter = 300
solver_options.initial_sto_reg_iter = 10 
ocp_solver = robotoc.OCPSolver(ocp=ocp, contact_sequence=contact_sequence, 
                               solver_options=solver_options, nthreads=4)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())
ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)

ocp_solver.mesh_refinement(t)

print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
ocp_solver.solve(t, q, v)
print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))
print(ocp_solver.get_solver_statistics())

# Display results
viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='gepetto')
viewer.set_contact_info(contact_frames, mu)
discretization = ocp_solver.get_time_discretization()
viewer.display(discretization.time_steps(), ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))