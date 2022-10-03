import robotoc
import numpy as np
import math


model_info = robotoc.RobotModelInfo()
model_info.urdf_path = '../anymal_b_simple_description/urdf/anymal.urdf'
model_info.base_joint_type = robotoc.BaseJointType.FloatingBase
baumgarte_time_step = 0.05
model_info.point_contacts = [robotoc.ContactModelInfo('LF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('LH_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RF_FOOT', baumgarte_time_step),
                             robotoc.ContactModelInfo('RH_FOOT', baumgarte_time_step)]
robot = robotoc.Robot(model_info)

dt = 0.02
jump_length = np.array([0.8, 0, 0])
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
q_ref[0:3] += jump_length
q_weight = np.array([1.0, 0., 0., 1.0, 1.0, 1.0, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001,
                     0.001, 0.001, 0.001])
v_weight = np.full(robot.dimv(), 1.0)
a_weight = np.full(robot.dimv(), 1.0e-06)
q_weight_impulse = np.array([0., 0., 0., 100., 100., 100., 
                      0.1, 0.1, 0.1, 
                      0.1, 0.1, 0.1,
                      0.1, 0.1, 0.1,
                      0.1, 0.1, 0.1])
v_weight_impulse = np.full(robot.dimv(), 1.0)
dv_weight_impulse = np.full(robot.dimv(), 1.0e-06)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_ref)
config_cost.set_q_weight(q_weight)
config_cost.set_q_weight_terminal(q_weight)
config_cost.set_q_weight_impulse(q_weight_impulse)
config_cost.set_v_weight(v_weight)
config_cost.set_v_weight_terminal(v_weight)
config_cost.set_v_weight_impulse(v_weight_impulse)
config_cost.set_dv_weight_impulse(dv_weight_impulse)
config_cost.set_a_weight(a_weight)
cost.push_back(config_cost)

# Create the constraints
constraints           = robotoc.Constraints(barrier_param=1.0e-03, fraction_to_boundary_rule=0.995)
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
friction_cone         = robotoc.FrictionCone(robot)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.push_back(friction_cone)


# Create the contact sequence
contact_sequence = robotoc.ContactSequence(robot)
mu = 0.7
friction_coefficients = {'LF_FOOT': mu, 'LH_FOOT': mu, 'RF_FOOT': mu, 'RH_FOOT': mu} 

robot.forward_kinematics(q_standing)
x3d0_LF = robot.frame_position('LF_FOOT')
x3d0_LH = robot.frame_position('LH_FOOT')
x3d0_RF = robot.frame_position('RF_FOOT')
x3d0_RH = robot.frame_position('RH_FOOT')
contact_positions = {'LF_FOOT': x3d0_LF, 'LH_FOOT': x3d0_LH, 'RF_FOOT': x3d0_RF, 'RH_FOOT': x3d0_RH} 

contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts(['LF_FOOT', 'LH_FOOT', 'RF_FOOT', 'RH_FOOT'])
contact_status_standing.set_contact_placements(contact_positions)
contact_status_standing.set_friction_coefficients(friction_coefficients)
contact_sequence.init(contact_status_standing)

contact_status_flying = robot.create_contact_status()
contact_sequence.push_back(contact_status_flying, t0+ground_time-0.3-0.05, sto=True)

contact_positions['LF_FOOT'] += jump_length
contact_positions['LH_FOOT'] += jump_length
contact_positions['RF_FOOT'] += jump_length
contact_positions['RH_FOOT'] += jump_length
contact_status_standing.set_contact_placements(contact_positions)
contact_sequence.push_back(contact_status_standing, t0+ground_time+flying_time+0.05, sto=True)

# you can check the contact sequence via 
# print(contact_sequence)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
sto_constraints = robotoc.STOConstraints(minimum_dwell_times=[0.15, 0.15, 0.65],
                                         barrier_param=1.0e-03, 
                                         fraction_to_boundary_rule=0.995)

T = t0 + flying_time + 2*ground_time
N = math.floor(T/dt) 
# Create the OCP with the STO problem
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  contact_sequence=contact_sequence, T=T, N=N)
# Create the OCP solver
solver_options = robotoc.SolverOptions()
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
solver_options.max_iter = 200
ocp_solver = robotoc.OCPSolver(ocp=ocp, solver_options=solver_options, nthreads=4)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

ocp_solver.init_constraints(t)
print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
ocp_solver.solve(t, q, v)
print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))
print(ocp_solver.get_solver_statistics())

# Plot results
kkt_data = ocp_solver.get_solver_statistics().kkt_error + [ocp_solver.KKT_error()] # append KKT after convergence
ts_data = ocp_solver.get_solver_statistics().ts + [contact_sequence.event_times()] # append ts after convergence

plot_ts = robotoc.utils.PlotConvergence()
plot_ts.ylim = [0., 1.5]
plot_ts.plot(kkt_data=kkt_data, ts_data=ts_data, fig_name='jump_sto', 
             save_dir='jump_sto_log')

plot_f = robotoc.utils.PlotContactForce(mu=mu)
plot_f.plot(f_data=ocp_solver.get_solution('f', 'WORLD'), 
            time_discretization=ocp_solver.get_time_discretization(), 
            fig_name='jump_sto_f', save_dir='jump_sto_log')

# Display results
viewer = robotoc.utils.TrajectoryViewer(model_info=model_info, viewer_type='gepetto')
viewer.set_contact_info(mu=mu)
viewer.display(ocp_solver.get_time_discretization(), 
               ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))