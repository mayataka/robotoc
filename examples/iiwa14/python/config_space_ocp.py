import robotoc
import numpy as np
import math


path_to_urdf = "../iiwa_description/urdf/iiwa14.urdf"
robot = robotoc.Robot(path_to_urdf)

# Change the limits from the default parameters.
robot.set_joint_effort_limit(np.full(robot.dimu(), 50))
robot.set_joint_velocity_limit(np.full(robot.dimv(), 0.5*math.pi))

# Create a cost function.
cost = robotoc.CostFunction()
config_cost = robotoc.ConfigurationSpaceCost(robot)
q_ref = np.array([0, 0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi, 0]) 
config_cost.set_q_ref(q_ref)
config_cost.set_q_weight(np.full(robot.dimv(), 10))
config_cost.set_qf_weight(np.full(robot.dimv(), 10))
config_cost.set_v_weight(np.full(robot.dimv(), 0.01))
config_cost.set_vf_weight(np.full(robot.dimv(), 0.01))
config_cost.set_a_weight(np.full(robot.dimv(), 0.01))
cost.push_back(config_cost)

# Create joint constraints.
constraints           = robotoc.Constraints()
joint_position_lower  = robotoc.JointPositionLowerLimit(robot)
joint_position_upper  = robotoc.JointPositionUpperLimit(robot)
joint_velocity_lower  = robotoc.JointVelocityLowerLimit(robot)
joint_velocity_upper  = robotoc.JointVelocityUpperLimit(robot)
joint_torques_lower   = robotoc.JointTorquesLowerLimit(robot)
joint_torques_upper   = robotoc.JointTorquesUpperLimit(robot)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)
constraints.set_barrier(1.0e-06)

# Create the OCP solver for unconstrained rigid-body systems.
T = 3.0
N = 60
nthreads = 4
t = 0.0
q = np.array([0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi]) 
v = np.zeros(robot.dimv())

ocp_solver = robotoc.UnconstrOCPSolver(robot, cost, constraints, T, N, nthreads)
ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)
ocp_solver.init_constraints()

print("----- Solves the OCP by Riccati recursion algorithm. -----")
num_iteration = 30
robotoc.utils.benchmark.convergence(ocp_solver, t, q, v, num_iteration)

# Solves the OCP by ParNMPC algorithm.
nthreads = 8
parnmpc_solver = robotoc.UnconstrParNMPCSolver(robot, cost, constraints, T, N, nthreads)
parnmpc_solver.set_solution("q", q)
parnmpc_solver.set_solution("v", v)
parnmpc_solver.init_constraints()

print("\n----- Solves the OCP by ParNMPC algorithm. -----")
num_iteration = 60
robotoc.utils.benchmark.convergence(parnmpc_solver, t, q, v, num_iteration)

viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, viewer_type='meshcat')
viewer.set_camera_transform_meshcat(camera_tf_vec=[0.5, -3.0, 0.0], zoom=2.0)
viewer.display((T/N), ocp_solver.get_solution('q'))