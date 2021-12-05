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
step_length = 0.15
step_height = 0.08
swing_time = 0.40
double_support_time = 0.05
t0 = 0.1
cycle = 1


# Create the cost function
cost = robotoc.CostFunction()
q_standing = np.array([0, 0, 0.4792, 0, 0, 0, 1, 
                       -0.1,  0.7, -1.0, 
                       -0.1, -0.7,  1.0, 
                        0.1,  0.7, -1.0, 
                        0.1, -0.7,  1.0])
q_ref = q_standing.copy()
q_ref[0] += step_length
q_weight = np.array([0, 10000, 10000, 10000, 10000, 10000, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001, 
                     0.001, 0.001, 0.001])
v_weight = np.full(robot.dimv(), 0.1)
v_weight[:6] = np.full(6, 100.0)
v_ref = np.full(robot.dimv(), 0.0)
v_ref[0] = 0.5 * step_length / (swing_time+double_support_time)
a_weight = np.full(robot.dimv(), 1.0e-03)
dvi_weight  = np.full(robot.dimv(), 1.0e-03)
config_cost = robotoc.ConfigurationSpaceCost(robot)
config_cost.set_q_ref(q_ref)
config_cost.set_v_ref(v_ref)
config_cost.set_q_weight(q_weight)
config_cost.set_qf_weight(q_weight)
config_cost.set_qi_weight(q_weight)
config_cost.set_v_weight(v_weight)
config_cost.set_vf_weight(v_weight)
config_cost.set_vi_weight(v_weight)
config_cost.set_a_weight(a_weight)
cost.push_back(config_cost)


LF_foot_ref = robotoc.TrottingSwingFootRef(0, 2, 1, step_length, step_height)
LH_foot_ref = robotoc.TrottingSwingFootRef(1, 3, 0, step_length, step_height)
RF_foot_ref = robotoc.TrottingSwingFootRef(2, 0, 3, step_length, step_height)
RH_foot_ref = robotoc.TrottingSwingFootRef(3, 1, 2, step_length, step_height)
LF_cost = robotoc.SwingFootCost(robot, 0, LF_foot_ref)
LH_cost = robotoc.SwingFootCost(robot, 1, LH_foot_ref)
RF_cost = robotoc.SwingFootCost(robot, 2, RF_foot_ref)
RH_cost = robotoc.SwingFootCost(robot, 3, RH_foot_ref)
point_ref_weight = np.array([0., 1.0e06, 1.0e04])
LF_cost.set_q_weight(point_ref_weight)
LH_cost.set_q_weight(point_ref_weight)
RF_cost.set_q_weight(point_ref_weight)
RH_cost.set_q_weight(point_ref_weight)
cost.push_back(LF_cost)
cost.push_back(LH_cost)
cost.push_back(RF_cost)
cost.push_back(RH_cost)

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
max_num_impulses = 2*cycle
contact_sequence = robotoc.ContactSequence(robot, max_num_impulses)

robot.forward_kinematics(q_standing)
q0_3d_LF = robot.frame_position(LF_foot_id)
q0_3d_LH = robot.frame_position(LH_foot_id)
q0_3d_RF = robot.frame_position(RF_foot_id)
q0_3d_RH = robot.frame_position(RH_foot_id)
contact_points = [q0_3d_LF, q0_3d_LH, q0_3d_RF, q0_3d_RH]

contact_status_standing = robot.create_contact_status()
contact_status_standing.activate_contacts([0, 1, 2, 3])
contact_status_standing.set_contact_points(contact_points)
contact_sequence.init_contact_sequence(contact_status_standing)

contact_status_lhrf_swing = robot.create_contact_status()
contact_status_lhrf_swing.activate_contacts([0, 3])
contact_status_lhrf_swing.set_contact_points(contact_points)
contact_sequence.push_back(contact_status=contact_status_lhrf_swing, 
                           switching_time=t0, sto=True)

contact_points[1][0] += 0.5 * step_length
contact_points[2][0] += 0.5 * step_length
contact_status_standing.set_contact_points(contact_points)
contact_sequence.push_back(contact_status_standing, t0+0.21, sto=True)

contact_status_lfrh_swing = robot.create_contact_status()
contact_status_lfrh_swing.activate_contacts([1, 2])
contact_status_lfrh_swing.set_contact_points(contact_points)
contact_sequence.push_back(contact_status_lfrh_swing, t0+0.42, sto=True)

contact_points[0][0] += step_length
contact_points[3][0] += step_length
contact_status_standing.set_contact_points(contact_points)
contact_sequence.push_back(contact_status_standing, t0+0.63, sto=True)

for i in range(cycle-1):
    t1 = t0 + (i+1)*0.84
    contact_status_lhrf_swing.set_contact_points(contact_points)
    contact_sequence.push_back(contact_status_lhrf_swing, t1, sto=True)

    contact_points[1][0] += step_length
    contact_points[2][0] += step_length
    contact_status_standing.set_contact_points(contact_points)
    contact_sequence.push_back(contact_status_standing, t1+0.21, sto=True)

    contact_status_lfrh_swing.set_contact_points(contact_points)
    contact_sequence.push_back(contact_status_lfrh_swing, 
                               t1+0.42, sto=True)

    contact_points[0][0] += step_length
    contact_points[3][0] += step_length
    contact_status_standing.set_contact_points(contact_points)
    contact_sequence.push_back(contact_status_standing, 
                               t1+0.63, sto=True)

# you can chech the contact sequence as 
# print(contact_sequence)

# Create the STO cost function. This is necessary even empty one to construct an OCP with a STO problem
sto_cost = robotoc.STOCostFunction()
# Create the STO constraints 
sto_constraints = robotoc.STOConstraints(2*max_num_impulses)
sto_constraints.set_minimum_dwell_times([0.02, 0.2, 0.02, 0.2, 0.02])
sto_constraints.set_barrier(1.0e-03)

T = t0 + cycle*(2*double_support_time+2*swing_time)
N = math.floor(T/dt) 
# Create the OCP with the STO problem
ocp = robotoc.OCP(robot=robot, cost=cost, constraints=constraints, 
                  sto_cost=sto_cost, sto_constraints=sto_constraints, 
                  T=T, N=N, max_num_each_discrete_events=max_num_impulses)
# Create the OCP solver
solver_options = robotoc.SolverOptions()
solver_options.print_level = 1
solver_options.kkt_tol_mesh = 0.1
solver_options.max_dt_mesh = T/N 
ocp_solver = robotoc.OCPSolver(ocp=ocp, contact_sequence=contact_sequence, 
                               solver_options=solver_options, nthreads=4)

# Initial time and intial state 
t = 0.
q = q_standing
v = np.zeros(robot.dimv())

ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
f_init = np.array([0.0, 0.0, 0.25*robot.total_weight()])
ocp_solver.set_solution("f", f_init)

ocp_solver.mesh_refinement(t)
# print("Initial KKT error: ", ocp_solver.KKT_error(t, q, v))
# ocp_solver.solve(t, q, v)
# print("KKT error after convergence: ", ocp_solver.KKT_error(t, q, v))

logger_kkt_ts = robotoc.utils.Logger(vars=['ts', 'KKT'], log_name='trotting_sto')
robotoc.utils.benchmark.convergence_sto(ocp_solver, t, q, v, num_iteration=50, 
                                        dt_tol_mesh=0.02, kkt_tol_mesh=0.1, 
                                        logger=logger_kkt_ts)

kkt_data = logger_kkt_ts.get_data('KKT')
ts_data = logger_kkt_ts.get_data('ts')
plot = robotoc.utils.PlotConvergence()
plot.ylim = [0., 1.2]
plot.plot(kkt_data=kkt_data, ts_data=ts_data, fig_name='trotting_sto', 
          save_dir=logger_kkt_ts.get_log_dir())

plot = robotoc.utils.PlotContactForce(mu=mu)
plot.plot(f_data=ocp_solver.get_solution('f', 'WORLD'), 
          t=ocp_solver.get_OCP_discretization().time_points(), 
          fig_name='trotting_sto_f', save_dir=logger_kkt_ts.get_log_dir())

viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='gepetto')
viewer.set_contact_info(contact_frames, mu)
ocp_discretization = ocp_solver.get_OCP_discretization()
viewer.display(ocp_discretization.time_steps(), ocp_solver.get_solution('q'), 
               ocp_solver.get_solution('f', 'WORLD'))