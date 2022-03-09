import robotoc_sim
import pybullet


class A1Simulator(robotoc_sim.LeggedSimulator):
    def __init__(self, path_to_urdf, time_step, start_time, end_time):
        super().__init__(path_to_urdf, time_step, start_time, end_time)

    @classmethod
    def get_state_from_pybullet(self, pybullet_robot, q, v):
        # Base
        basePos, baseOrn = pybullet.getBasePositionAndOrientation(pybullet_robot)
        q[0:3] = basePos
        q[3:7] = baseOrn
        # FL
        q[7] = pybullet.getJointState(pybullet_robot, 7)[0]
        q[8] = pybullet.getJointState(pybullet_robot, 9)[0]
        q[9] = pybullet.getJointState(pybullet_robot, 10)[0]
        # FR
        q[10] = pybullet.getJointState(pybullet_robot, 2)[0]
        q[11] = pybullet.getJointState(pybullet_robot, 4)[0]
        q[12] = pybullet.getJointState(pybullet_robot, 5)[0]
        # RF
        q[13] = pybullet.getJointState(pybullet_robot, 17)[0]
        q[14] = pybullet.getJointState(pybullet_robot, 19)[0]
        q[15] = pybullet.getJointState(pybullet_robot, 20)[0]
        # RH
        q[16] = pybullet.getJointState(pybullet_robot, 12)[0]
        q[17] = pybullet.getJointState(pybullet_robot, 14)[0]
        q[18] = pybullet.getJointState(pybullet_robot, 15)[0]

        # Base
        baseVel, baseAngVel = self.get_body_local_velocity(pybullet_robot)
        v[0:3] = baseVel
        v[3:6] = baseAngVel
        # FL
        v[6] = pybullet.getJointState(pybullet_robot, 7)[1]
        v[7] = pybullet.getJointState(pybullet_robot, 9)[1]
        v[8] = pybullet.getJointState(pybullet_robot, 10)[1]
        # FR
        v[9] = pybullet.getJointState(pybullet_robot, 2)[1]
        v[10] = pybullet.getJointState(pybullet_robot, 4)[1]
        v[11] = pybullet.getJointState(pybullet_robot, 5)[1]
        # RF
        v[12] = pybullet.getJointState(pybullet_robot, 17)[1]
        v[13] = pybullet.getJointState(pybullet_robot, 19)[1]
        v[14] = pybullet.getJointState(pybullet_robot, 20)[1]
        # RH
        v[15] = pybullet.getJointState(pybullet_robot, 12)[1]
        v[16] = pybullet.getJointState(pybullet_robot, 14)[1]
        v[17] = pybullet.getJointState(pybullet_robot, 15)[1]

    @classmethod
    def apply_control_input_to_pybullet(self, pybullet_robot, u):
        mode = pybullet.TORQUE_CONTROL
        # LFL
        pybullet.setJointMotorControl2(pybullet_robot, 7, controlMode=mode, force=u[0])
        pybullet.setJointMotorControl2(pybullet_robot, 9, controlMode=mode, force=u[1])
        pybullet.setJointMotorControl2(pybullet_robot, 10, controlMode=mode, force=u[2])
        # FR 
        pybullet.setJointMotorControl2(pybullet_robot, 2, controlMode=mode, force=u[3])
        pybullet.setJointMotorControl2(pybullet_robot, 4, controlMode=mode, force=u[4])
        pybullet.setJointMotorControl2(pybullet_robot, 5, controlMode=mode, force=u[5])
        # RF
        pybullet.setJointMotorControl2(pybullet_robot, 17, controlMode=mode, force=u[6])
        pybullet.setJointMotorControl2(pybullet_robot, 19, controlMode=mode, force=u[7])
        pybullet.setJointMotorControl2(pybullet_robot, 20, controlMode=mode, force=u[8])
        # RH
        pybullet.setJointMotorControl2(pybullet_robot, 12, controlMode=mode, force=u[9])
        pybullet.setJointMotorControl2(pybullet_robot, 14, controlMode=mode, force=u[10])
        pybullet.setJointMotorControl2(pybullet_robot, 15, controlMode=mode, force=u[11])


# pin: FL_hip, FL_thigh, FL_calf, FR_hip, FR_thigh, FR_calf (6, 7, 8, 9, 10, 11)
# pin: RL_hip, RL_thigh, RL_calf, RR_hip, RR_thigh, RR_calf (12, 13, 14, 15, 16, 17)

    @classmethod
    def init_pybullet_robot(self, pybullet_robot, q):
        # Base
        pybullet.resetBasePositionAndOrientation(pybullet_robot, q[0:3], q[3:7])
        # FR
        pybullet.resetJointState(pybullet_robot, 2, q[10])
        pybullet.resetJointState(pybullet_robot, 4, q[11])
        pybullet.resetJointState(pybullet_robot, 5, q[12])
        # FL
        pybullet.resetJointState(pybullet_robot, 7, q[7])
        pybullet.resetJointState(pybullet_robot, 9, q[8])
        pybullet.resetJointState(pybullet_robot, 10, q[9])
        # RR 
        pybullet.resetJointState(pybullet_robot, 12, q[16])
        pybullet.resetJointState(pybullet_robot, 14, q[17])
        pybullet.resetJointState(pybullet_robot, 15, q[18])
        # RL 
        pybullet.resetJointState(pybullet_robot, 17, q[13])
        pybullet.resetJointState(pybullet_robot, 19, q[14])
        pybullet.resetJointState(pybullet_robot, 20, q[15])
        # Turn off the velocity and position control effects
        joints = [1, 2, 4, 5, 7, 9, 10, 12, 14, 15, 17, 19, 20]
        for j in joints:
            pybullet.setJointMotorControl2(pybullet_robot, j, 
                                           controlMode=pybullet.VELOCITY_CONTROL, 
                                           force=0.)
            pybullet.setJointMotorControl2(pybullet_robot, j, 
                                           controlMode=pybullet.POSITION_CONTROL, 
                                           force=0.)