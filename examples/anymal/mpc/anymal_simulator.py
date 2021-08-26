import idocp_sim
import pybullet


class ANYmalSimulator(idocp_sim.QuadrupedalSimulator):
    def __init__(self, path_to_urdf, time_step, start_time, end_time):
        super().__init__(path_to_urdf, time_step, start_time, end_time)

    @classmethod
    def get_state_from_pybullet(self, pybullet_robot, q, v):
        # Base
        basePos, baseOrn = pybullet.getBasePositionAndOrientation(pybullet_robot)
        q[0:3] = basePos
        q[3:7] = baseOrn
        # LF 
        q[7] = pybullet.getJointState(pybullet_robot, 1)[0]
        q[8] = pybullet.getJointState(pybullet_robot, 2)[0]
        q[9] = pybullet.getJointState(pybullet_robot, 3)[0]
        # LH 
        q[10] = pybullet.getJointState(pybullet_robot, 11)[0]
        q[11] = pybullet.getJointState(pybullet_robot, 12)[0]
        q[12] = pybullet.getJointState(pybullet_robot, 13)[0]
        # RF
        q[13] = pybullet.getJointState(pybullet_robot, 6)[0]
        q[14] = pybullet.getJointState(pybullet_robot, 7)[0]
        q[15] = pybullet.getJointState(pybullet_robot, 8)[0]
        # RH
        q[16] = pybullet.getJointState(pybullet_robot, 16)[0]
        q[17] = pybullet.getJointState(pybullet_robot, 17)[0]
        q[18] = pybullet.getJointState(pybullet_robot, 18)[0]
        # Base
        baseVel, baseAngVel = pybullet.getBaseVelocity(pybullet_robot)
        v[0:3] = baseVel
        v[3:6] = baseAngVel
        # LF 
        v[6] = pybullet.getJointState(pybullet_robot, 1)[1]
        v[7] = pybullet.getJointState(pybullet_robot, 2)[1]
        v[8] = pybullet.getJointState(pybullet_robot, 3)[1]
        # LH 
        v[9] = pybullet.getJointState(pybullet_robot, 11)[1]
        v[10] = pybullet.getJointState(pybullet_robot, 12)[1]
        v[11] = pybullet.getJointState(pybullet_robot, 13)[1]
        # RF
        v[12] = pybullet.getJointState(pybullet_robot, 6)[1]
        v[13] = pybullet.getJointState(pybullet_robot, 7)[1]
        v[14] = pybullet.getJointState(pybullet_robot, 8)[1]
        # RH
        v[15] = pybullet.getJointState(pybullet_robot, 16)[1]
        v[16] = pybullet.getJointState(pybullet_robot, 17)[1]
        v[17] = pybullet.getJointState(pybullet_robot, 18)[1]

    @classmethod
    def apply_control_input_to_pybullet(self, pybullet_robot, u):
        mode = pybullet.TORQUE_CONTROL
        # LF 
        pybullet.setJointMotorControl2(pybullet_robot, 1, controlMode=mode, force=u[0])
        pybullet.setJointMotorControl2(pybullet_robot, 2, controlMode=mode, force=u[1])
        pybullet.setJointMotorControl2(pybullet_robot, 3, controlMode=mode, force=u[2])
        # LH 
        pybullet.setJointMotorControl2(pybullet_robot, 11, controlMode=mode, force=u[3])
        pybullet.setJointMotorControl2(pybullet_robot, 12, controlMode=mode, force=u[4])
        pybullet.setJointMotorControl2(pybullet_robot, 13, controlMode=mode, force=u[5])
        # RF
        pybullet.setJointMotorControl2(pybullet_robot, 6, controlMode=mode, force=u[6])
        pybullet.setJointMotorControl2(pybullet_robot, 7, controlMode=mode, force=u[7])
        pybullet.setJointMotorControl2(pybullet_robot, 8, controlMode=mode, force=u[8])
        # RH
        pybullet.setJointMotorControl2(pybullet_robot, 16, controlMode=mode, force=u[9])
        pybullet.setJointMotorControl2(pybullet_robot, 17, controlMode=mode, force=u[10])
        pybullet.setJointMotorControl2(pybullet_robot, 18, controlMode=mode, force=u[11])

    @classmethod
    def init_pybullet_robot(self, pybullet_robot, q):
        # Base
        pybullet.resetBasePositionAndOrientation(pybullet_robot, q[0:3], q[3:7])
        # LF 
        pybullet.resetJointState(pybullet_robot, 1, q[7])
        pybullet.resetJointState(pybullet_robot, 2, q[8])
        pybullet.resetJointState(pybullet_robot, 3, q[9])
        # LH 
        pybullet.resetJointState(pybullet_robot, 11, q[10])
        pybullet.resetJointState(pybullet_robot, 12, q[11])
        pybullet.resetJointState(pybullet_robot, 13, q[12])
        # RF 
        pybullet.resetJointState(pybullet_robot, 6, q[13])
        pybullet.resetJointState(pybullet_robot, 7, q[14])
        pybullet.resetJointState(pybullet_robot, 8, q[15])
        # RH 
        pybullet.resetJointState(pybullet_robot, 16, q[16])
        pybullet.resetJointState(pybullet_robot, 17, q[17])
        pybullet.resetJointState(pybullet_robot, 18, q[18])
        # Turn off the velocity and position control effects
        joints = [1, 2, 3, 11, 12, 13, 6, 7, 8, 16, 17, 18]
        for j in joints:
            pybullet.setJointMotorControl2(pybullet_robot, j, 
                                           controlMode=pybullet.VELOCITY_CONTROL, 
                                           force=0.)
            pybullet.setJointMotorControl2(pybullet_robot, j, 
                                           controlMode=pybullet.POSITION_CONTROL, 
                                           force=0.)