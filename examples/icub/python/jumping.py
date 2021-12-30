import robotoc
import numpy as np
import math


L_foot_id = 26
R_foot_id = 46
contact_frames = [L_foot_id, R_foot_id]
contact_types = [robotoc.ContactType.SurfaceContact, robotoc.ContactType.SurfaceContact]
path_to_urdf = '../icub_description/urdf/icub.urdf'
baumgarte_time_step = 0.04
robot = robotoc.Robot(path_to_urdf, robotoc.BaseJointType.FloatingBase, 
                      contact_frames, contact_types, baumgarte_time_step)
print(robot)
q = np.zeros(robot.dimq())
q[6] = 1.0

# cost = robotoc.CostFunction()
# q_standing = np.array([0, 0, 0.4792, 0, 0, 0, 1, 
#                        -0.1,  0.7, -1.0, 
#                        -0.1, -0.7,  1.0, 
#                         0.1,  0.7, -1.0, 
#                         0.1, -0.7,  1.0])

#         <joint name="root_joint" value="0. 0. 0.60 0. 0. 1. 0." />
# robot.forward_kinematics(q)


viewer = robotoc.utils.TrajectoryViewer(path_to_urdf=path_to_urdf, 
                                        base_joint_type=robotoc.BaseJointType.FloatingBase,
                                        viewer_type='meshcat')
viewer.display(0.01, [q for i in range(50)])