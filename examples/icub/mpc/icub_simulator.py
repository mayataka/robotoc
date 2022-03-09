import robotoc_sim
import pybullet


class iCubSimulator(robotoc_sim.LeggedSimulator):
    def __init__(self, path_to_urdf, time_step, start_time, end_time):
        super().__init__(path_to_urdf, time_step, start_time, end_time)

    @classmethod
    def get_state_from_pybullet(self, pybullet_robot, q, v):
        # Base
        basePos, baseOrn = pybullet.getBasePositionAndOrientation(pybullet_robot)
        q[0:3] = basePos
        q[3:7] = baseOrn
        # L
        q[7]  = pybullet.getJointState(pybullet_robot, 1)[0]
        q[8]  = pybullet.getJointState(pybullet_robot, 2)[0]
        q[9]  = pybullet.getJointState(pybullet_robot, 4)[0]
        q[10] = pybullet.getJointState(pybullet_robot, 5)[0]
        q[11] = pybullet.getJointState(pybullet_robot, 6)[0]
        q[12] = pybullet.getJointState(pybullet_robot, 7)[0]
        # R
        q[13] = pybullet.getJointState(pybullet_robot, 12)[0]
        q[14] = pybullet.getJointState(pybullet_robot, 13)[0]
        q[15] = pybullet.getJointState(pybullet_robot, 15)[0]
        q[16] = pybullet.getJointState(pybullet_robot, 16)[0]
        q[17] = pybullet.getJointState(pybullet_robot, 17)[0]
        q[18] = pybullet.getJointState(pybullet_robot, 18)[0]

        # Base
        baseVel, baseAngVel = self.get_body_local_velocity(pybullet_robot)
        v[0:3] = baseVel
        v[3:6] = baseAngVel
        # L
        v[6]  = pybullet.getJointState(pybullet_robot, 1)[1]
        v[7]  = pybullet.getJointState(pybullet_robot, 2)[1]
        v[8]  = pybullet.getJointState(pybullet_robot, 4)[1]
        v[9]  = pybullet.getJointState(pybullet_robot, 5)[1]
        v[10] = pybullet.getJointState(pybullet_robot, 6)[1]
        v[11] = pybullet.getJointState(pybullet_robot, 7)[1]
        # R
        v[12] = pybullet.getJointState(pybullet_robot, 12)[1]
        v[13] = pybullet.getJointState(pybullet_robot, 13)[1]
        v[14] = pybullet.getJointState(pybullet_robot, 15)[1]
        v[15] = pybullet.getJointState(pybullet_robot, 16)[1]
        v[16] = pybullet.getJointState(pybullet_robot, 17)[1]
        v[17] = pybullet.getJointState(pybullet_robot, 18)[1]

    @classmethod
    def apply_control_input_to_pybullet(self, pybullet_robot, u):
        mode = pybullet.TORQUE_CONTROL
        # L
        pybullet.setJointMotorControl2(pybullet_robot, 1, controlMode=mode, force=u[0])
        pybullet.setJointMotorControl2(pybullet_robot, 2, controlMode=mode, force=u[1])
        pybullet.setJointMotorControl2(pybullet_robot, 4, controlMode=mode, force=u[2])
        pybullet.setJointMotorControl2(pybullet_robot, 5, controlMode=mode, force=u[3])
        pybullet.setJointMotorControl2(pybullet_robot, 6, controlMode=mode, force=u[4])
        pybullet.setJointMotorControl2(pybullet_robot, 7, controlMode=mode, force=u[5])
        # R
        pybullet.setJointMotorControl2(pybullet_robot, 12, controlMode=mode, force=u[6])
        pybullet.setJointMotorControl2(pybullet_robot, 13, controlMode=mode, force=u[7])
        pybullet.setJointMotorControl2(pybullet_robot, 15, controlMode=mode, force=u[8])
        pybullet.setJointMotorControl2(pybullet_robot, 16, controlMode=mode, force=u[9])
        pybullet.setJointMotorControl2(pybullet_robot, 17, controlMode=mode, force=u[10])
        pybullet.setJointMotorControl2(pybullet_robot, 18, controlMode=mode, force=u[11])

    @classmethod
    def init_pybullet_robot(self, pybullet_robot, q):
        # Base
        pybullet.resetBasePositionAndOrientation(pybullet_robot, q[0:3], q[3:7])
        # L
        pybullet.resetJointState(pybullet_robot, 1, q[7])
        pybullet.resetJointState(pybullet_robot, 2, q[8])
        pybullet.resetJointState(pybullet_robot, 4, q[9])
        pybullet.resetJointState(pybullet_robot, 5, q[10])
        pybullet.resetJointState(pybullet_robot, 6, q[11])
        pybullet.resetJointState(pybullet_robot, 7, q[12])
        # R
        pybullet.resetJointState(pybullet_robot, 12, q[13])
        pybullet.resetJointState(pybullet_robot, 13, q[14])
        pybullet.resetJointState(pybullet_robot, 15, q[15])
        pybullet.resetJointState(pybullet_robot, 16, q[16])
        pybullet.resetJointState(pybullet_robot, 17, q[17])
        pybullet.resetJointState(pybullet_robot, 18, q[18])

        for i in [1, 2, 4, 5, 6, 7, 12, 13, 15, 16, 17, 18]:
            pybullet.setJointMotorControl2(pybullet_robot, i, 
                                           controlMode=pybullet.VELOCITY_CONTROL, 
                                           force=0.)
            pybullet.setJointMotorControl2(pybullet_robot, i, 
                                           controlMode=pybullet.POSITION_CONTROL, 
                                           force=0.)
