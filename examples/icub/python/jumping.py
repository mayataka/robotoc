import robotoc
import numpy as np
import math


L_foot_id = 26
R_foot_id = 46
contact_frames = [L_foot_id, R_foot_id]
contact_types = [robotoc.ContactType.SurfaceContact, robotoc.ContactType.SurfaceContact]
path_to_urdf = '../icub_description/urdf/icub.urdf'
baumgarte_time_step = 0.04
contact_inv_damping = 1.0e-12
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step, contact_inv_damping)

dt = 0.02
jump_length = 0.5
flying_up_time = 0.15
flying_down_time = flying_up_time
flying_time = flying_up_time + flying_down_time
ground_time = 0.7
t0 = 0.

cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.592, 0, 0, 1, 0,
                       0.20944, 0.08727, 0, -0.1745, -0.0279, -0.08726, # right leg
                       0.20944, 0.08727, 0, -0.1745, -0.0279, -0.08726, # left leg
                       0, 0, 0, # torso
                       0, 0.35, 0.5, 0.5, 0, 0, 0, # left arm 
                       0, 0.35, 0.5, 0.5, 0, 0, 0]) # right arm 
q_weight = np.array([0, 0, 0, 250000, 250000, 250000, 
                     0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 
                     0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 
                     0.0001, 0.0001, 0.0001,
                     0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
                     0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001])
v_weight = np.array([100, 100, 100, 100, 100, 100, 
                     1, 1, 1, 1, 1, 1, 
                     1, 1, 1, 1, 1, 1, 
                     1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 1, 1, 1])
u_weight = np.full(robot.dimu(), 1.0e-01)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_standing)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(q_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(v_weight)
config_cost.set_u_weight(u_weight)
cost.push_back(config_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier=1.0e-03, fraction_to_boundary_rule=0.995)
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

# Create the contact sequence
max_num_impulses = 3
contact_sequence = robotoc.ContactSequence(robot, max_num_impulses)

robot.forward_kinematics(q_standing)
x3d0_L = robot.frame_placement(L_foot_id)
x3d0_R = robot.frame_placement(R_foot_id)
contact_placements = [x3d0_L, x3d0_R]

contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts([0, 1])
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_flying = robot.create_contact_status()
contact_sequence.push_back(contact_status_flying, t0+ground_time-0.3, sto=True)

contact_placements[0].trans[0] += jump_length
contact_placements[1].trans[0] += jump_length
contact_status_standing.set_contact_placements(contact_placements)
contact_sequence.push_back(contact_status_standing, t0+ground_time+flying_time-0.1, sto=True)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
sto_constraints = robotoc.STOConstraints(max_num_switches=2*max_num_impulses, 
                                         min_dt=[0.15, 0.15, 0.65],
                                         barrier=1.0e-03, 
                                         fraction_to_boundary_rule=0.995)

T = t0 + flying_time + 2*ground_time
N = math.floor(T/dt) 
# Create the OCP with the STO problem
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  T=T, N=N, max_num_each_discrete_events=max_num_impulses)
# Create the OCP solver
solver_options = robotoc.SolverOptions()
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
solver_options.max_iter = 200
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

viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='gepetto')
discretization = ocp_solver.get_time_discretization()
viewer.display(discretization.time_steps(), ocp_solver.get_solution('q'))