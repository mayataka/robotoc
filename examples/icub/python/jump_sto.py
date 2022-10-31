import robotoc
import numpy as np
import math


model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../icub_description/urdf/icub.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.surface_contacts = [robotoc.ContactModelInfo('l_sole', baumgarte_time_step),
                               robotoc.ContactModelInfo('r_sole', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

dt = 0.02
jump_length = np.array([0.5, 0, 0])
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
q_ref[0:3] += jump_length
q_weight = np.array([0, 1, 1, 100, 100, 100, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 1, 1,
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])
q_weight_terminal = q_weight
v_weight = np.full(robot.dimv(), 1.0e-03)
a_weight = np.full(robot.dimv(), 1.0e-05)
q_weight_impact = 1.0 * q_weight
v_weight_impact = 1.0 * v_weight
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_q_weight_terminal(q_weight)
config_cost.set_q_weight_impact(q_weight_impact)
config_cost.set_v_weight(v_weight)
config_cost.set_v_weight_terminal(v_weight)
config_cost.set_v_weight_impact(v_weight_impact)
config_cost.set_a_weight(a_weight)
cost.add("config_cost", config_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier_param=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
# contact_wrench_cone  = robotoc.ContactWrenchCone(robot, 0.1, 0.05)
friction_cone  = robotoc.FrictionCone(robot)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
# constraints.push_back(contact_wrench_cone)
constraints.push_back(friction_cone)

# Create the contact sequence
contact_sequence = robotoc.ContactSequence(robot)
mu = 0.6
friction_coefficients = {'l_sole': mu, 'r_sole': mu} 

robot.forward_kinematics(q_standing)
x3d0_L = robot.frame_placement('l_sole')
x3d0_R = robot.frame_placement('r_sole')
contact_placements = {'l_sole': x3d0_L, 'r_sole': x3d0_R} 

contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts(['l_sole', 'r_sole'])
contact_status_standing.set_contact_placements(contact_placements)
contact_status_standing.set_friction_coefficients(friction_coefficients)
contact_sequence.init(contact_status_standing)

contact_status_flying = robot.create_contact_status()
contact_sequence.push_back(contact_status_flying, t0+ground_time, sto=True)

contact_placements['l_sole'].trans = contact_placements['l_sole'].trans + jump_length
contact_placements['r_sole'].trans = contact_placements['r_sole'].trans + jump_length 
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.push_back(contact_status_standing, t0+ground_time+flying_time, sto=True)

contact_sequence.push_back(contact_status_flying, t0+2*ground_time+flying_time, sto=True)

contact_placements['l_sole'].trans = contact_placements['l_sole'].trans + jump_length 
contact_placements['r_sole'].trans = contact_placements['r_sole'].trans + jump_length 
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.push_back(contact_status_standing, t0+2*ground_time+2*flying_time, sto=True)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
sto_constraints = robotoc.STOConstraints(minimum_dwell_times=[0.6, 0.2, 0.6, 0.2, 0.6],
                                         barrier_param=1.0e-03, 
                                         fraction_to_boundary_rule=0.995)

T = t0 + 2*flying_time + 3*ground_time
N = math.floor(T/dt) 
# Create the OCP with the STO problem
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  contact_sequence=contact_sequence, T=T, N=N)
# Create the OCP solver
solver_options = robotoc.SolverOptions()
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
solver_options.max_iter = 350
solver_options.initial_sto_reg_iter = 10 
solver_options.nthreads = 4
ocp_solver = robotoc.OCPSolver(ocp=ocp, solver_options=solver_options)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())
ocp_solver.discretize(t)
ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)

ocp_solver.init_constraints()
print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
ocp_solver.solve(t, q, v)
print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))
print(ocp_solver.get_solver_statistics())

# Display results
viewer = robotoc.utils.TrajectoryViewer(model_info=model_info, viewer_type='gepetto')
viewer.set_contact_info(mu=mu)
viewer.display(ocp_solver.get_time_discretization(), 
               ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))