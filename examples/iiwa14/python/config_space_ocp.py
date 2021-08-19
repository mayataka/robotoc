import idocp
import numpy as np
import math


path_to_urdf = "../iiwa_description/urdf/iiwa14.urdf"
robot = idocp.Robot(path_to_urdf)

# Change the limits from the default parameters.
robot.set_joint_effort_limit(np.full(robot.dimu(), 50))
robot.set_joint_velocity_limit(np.full(robot.dimv(), 0.5*math.pi))

# Create a cost function.
cost = idocp.CostFunction()
config_cost = idocp.ConfigurationSpaceCost(robot)
q_ref = np.array([0, 0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi, 0]) 
config_cost.set_q_ref(q_ref)
config_cost.set_q_weight(np.full(robot.dimv(), 10))
config_cost.set_qf_weight(np.full(robot.dimv(), 10))
config_cost.set_v_weight(np.full(robot.dimv(), 0.01))
config_cost.set_vf_weight(np.full(robot.dimv(), 0.01))
config_cost.set_a_weight(np.full(robot.dimv(), 0.01))
cost.push_back(config_cost)

# Create joint constraints.
constraints           = idocp.Constraints()
joint_position_lower  = idocp.JointPositionLowerLimit(robot)
joint_position_upper  = idocp.JointPositionUpperLimit(robot)
joint_velocity_lower  = idocp.JointVelocityLowerLimit(robot)
joint_velocity_upper  = idocp.JointVelocityUpperLimit(robot)
joint_torques_lower   = idocp.JointTorquesLowerLimit(robot)
joint_torques_upper   = idocp.JointTorquesUpperLimit(robot)
constraints.push_back(joint_position_lower)
constraints.push_back(joint_position_upper)
constraints.push_back(joint_velocity_lower)
constraints.push_back(joint_velocity_upper)
constraints.push_back(joint_torques_lower)
constraints.push_back(joint_torques_upper)

# Create the OCP solver for unconstrained rigid-body systems.
T = 3.0
N = 60
nthreads = 4
t = 0.0
q = np.array([0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi, 0, 0.5*math.pi]) 
v = np.zeros(robot.dimv())
ocp_solver = idocp.UnconstrOCPSolver(robot, cost, constraints, T, N, nthreads)
# ocp_solver = idocp.UnconstrParNMPCSolver(robot, cost, constraints, T, N, nthreads)

# Solves the OCP.
ocp_solver.set_solution("q", q)
ocp_solver.set_solution("v", v)

ocp_solver.init_constraints()

num_iteration = 30
idocp.benchmark_convergence(ocp_solver, t, q, v, num_iteration)