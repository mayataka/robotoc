import pybullet
import time
import abc
import numpy as np
import os


class CameraSettings(object):
    def __init__(self, camera_distance, camera_yaw, camera_pitch, camera_target_pos):
        self.camera_distance = camera_distance
        self.camera_yaw = camera_yaw
        self.camera_pitch = camera_pitch
        self.camera_target_pos = camera_target_pos


class TerrainSettings(object):
    def __init__(self, height_range: float=1.0, num_heightfield_rows: int=100, num_heightfield_columns: int=100, from_urdf=True):
        self.height_range = height_range
        self.num_heightfield_rows = num_heightfield_rows
        self.num_heightfield_columns = num_heightfield_columns
        self.terrain_urdf = None
        if from_urdf:
            self.terrain_urdf = os.path.join(os.path.dirname(__file__), "rsc/terrain.urdf")

    def create_random_terrain_shape(self):
        heightfield_data = np.zeros(shape=[self.num_heightfield_columns, self.num_heightfield_rows], dtype=np.float)
        for i in range(int(self.num_heightfield_columns/2)):
            for j in range(int(self.num_heightfield_rows)):
                n1 = 0
                n2 = 0
                if j > 0:
                    n1 = heightfield_data[i, j-1]
                if i > 0:
                    n2 = heightfield_data[i-1, j]
                else:
                    n2 = n1
                noise = np.random.uniform(-1.0, 1.0)
                heightfield_data[i, j] = (n1+n2)/2 + noise
        max = np.max(heightfield_data)
        min = np.min(heightfield_data)
        if max - min > 0:
            heightfield_data = (2.0 * self.height_range / (max-min)) * heightfield_data 
        heightfield_data_inv = heightfield_data[::-1,:]
        heightfield_data_2 = np.concatenate((heightfield_data_inv, heightfield_data))
        col, row = heightfield_data_2.shape
        heightfield_data_2 = heightfield_data_2.reshape(-1)
        terrain_shape = pybullet.createCollisionShape(shapeType=pybullet.GEOM_HEIGHTFIELD, 
                                                      heightfieldData=heightfield_data_2, 
                                                      meshScale=[0.5, 0.5, 1.0], 
                                                      numHeightfieldRows=row, 
                                                      numHeightfieldColumns=col)
        return terrain_shape


class ContactInfo(object):
    def __init__(self, link_id, link_name, contact_normal, contact_force):
        self.link_id = link_id
        self.link_name = link_name
        self.contact_normal = contact_normal
        self.contact_force = contact_force


class LeggedSimulator(metaclass=abc.ABCMeta):
    def __init__(self, urdf_path, time_step: float):
        self.urdf_path = urdf_path
        self.time_step = time_step
        self.camera_settings = None
        self.terrain_settings = None
        self.robot = None
        self.world = None
        self.t = 0
        self.q = None 
        self.v = None 
        self.u = None 
        self.base_lin_vel_world_prev = np.zeros(3)
    
    def set_camera_settings(self, camera_settings: CameraSettings):
        self.camera_settings = camera_settings

    def set_terrain_settings(self, terrain_settings: TerrainSettings):
        self.terrain_settings = terrain_settings

    def get_time(self):
        return self.t

    @abc.abstractmethod
    def get_joint_id_map(self):
        return NotImplementedError()

    def get_base_rotation_matrix(self):
        base_pos, base_orn = pybullet.getBasePositionAndOrientation(self.robot)
        R = np.reshape(pybullet.getMatrixFromQuaternion(base_orn), [3, 3]) 
        return R

    def get_body_world_velocity(self):
        base_lin_vel_world, base_ang_vel_world = pybullet.getBaseVelocity(self.robot)
        return np.array(base_lin_vel_world), np.array(base_ang_vel_world)

    def get_body_local_velocity(self):
        R = self.get_base_rotation_matrix()
        base_lin_vel_world, base_ang_vel_world = self.get_body_world_velocity()
        base_lin_vel_local = R.T @ np.array(base_lin_vel_world)
        base_ang_vel_local = R.T @ np.array(base_ang_vel_world)
        return base_lin_vel_local, base_ang_vel_local

    def get_imu_measurements(self):
        _, base_ang_vel_local = self.get_body_local_velocity()
        base_lin_vel_world, _ = self.get_body_world_velocity()
        base_lin_acc_world = (base_lin_vel_world - self.base_lin_vel_world_prev) / self.time_step + np.array([0, 0, 9.81])
        R = self.get_base_rotation_matrix()
        base_lin_acc_local = R.T @ base_lin_acc_world
        return base_ang_vel_local, base_lin_acc_local 

    def get_joint_measurements(self):
        q, v = self.get_state(self.robot)
        num_joints = len(self.get_joint_id_map())
        return q[-num_joints:], v[-num_joints:], self.u

    def get_state(self):
        joint_id_map = self.get_joint_id_map()
        num_joints = len(joint_id_map)
        q = np.zeros(7+num_joints)
        base_pos, base_orn = pybullet.getBasePositionAndOrientation(self.robot)
        q[0:3] = base_pos
        q[3:7] = base_orn
        for i in range(num_joints):
            q[7+i] = pybullet.getJointState(self.robot, joint_id_map[i])[0]
        v = np.zeros(6+num_joints)
        base_lin_vel, base_ang_vel = self.get_body_local_velocity()
        v[0:3] = base_lin_vel
        v[3:6] = base_ang_vel
        for i in range(num_joints):
            v[6+i] = pybullet.getJointState(self.robot, joint_id_map[i])[1]
        return q, v

    def get_contact_info(self):
        contact_points = pybullet.getContactPoints(self.robot, self.plane)
        contact_info = []
        for e in contact_points:
            link_id = e[3]
            link_name = pybullet.getJointInfo(self.robot, link_id)[1].decode()
            contact_normal = np.array(e[7])
            contact_force = e[9] * contact_normal
            contact_info.append(ContactInfo(link_id, link_name, contact_normal, contact_force))
        return contact_info

    def apply_control_input(self, u: np.ndarray):
        joint_id_map = self.get_joint_id_map()
        num_joints = len(joint_id_map)
        assert u.shape[0] == num_joints
        mode = pybullet.TORQUE_CONTROL
        for i in range(num_joints):
            pybullet.setJointMotorControl2(self.robot, joint_id_map[i], controlMode=mode, force=u[i])

    def init_state(self, q: np.ndarray):
        joint_id_map = self.get_joint_id_map()
        num_joints = len(joint_id_map)
        assert q.shape[0] == 7 + num_joints
        pybullet.resetBasePositionAndOrientation(self.robot, q[0:3], q[3:7])
        for i in range(num_joints):
            pybullet.resetJointState(self.robot, joint_id_map[i], q[7+i])
        for j in joint_id_map:
            pybullet.setJointMotorControl2(self.robot, j, 
                                           controlMode=pybullet.VELOCITY_CONTROL, 
                                           force=0.)
            pybullet.setJointMotorControl2(self.robot, j, 
                                           controlMode=pybullet.POSITION_CONTROL, 
                                           force=0.)

    def init_simulation(self, t0: float, q0: np.ndarray):
        pybullet.connect(pybullet.GUI)
        pybullet.setGravity(0, 0, -9.81)
        pybullet.setTimeStep(self.time_step)
        if self.terrain_settings is not None:
            if self.terrain_settings.terrain_urdf is not None:
                self.world = pybullet.loadURDF(self.terrain_settings.terrain_urdf)
            else:
                terrain_shape = self.terrain_settings.create_random_terrain_shape()
                self.world = pybullet.createMultiBody(0, terrain_shape)
                pybullet.changeVisualShape(self.world, -1, rgbaColor=[1.0, 1.0, 1.0, 1.0])
                pybullet.resetBasePositionAndOrientation(self.world, [0., 0., 0.], [0., 0., 0., 1.])
                ray_test = pybullet.rayTest([q0[0], q0[1], self.height_range], 
                                            [q0[0], q0[1], -self.height_range])
                ray_hit_position = ray_test[0][3]
                pybullet.resetBasePositionAndOrientation(self.world, [0., 0., -ray_hit_position[2]], [0., 0., 0., 1.])
        else:
            import pybullet_data
            pybullet.setAdditionalSearchPath(pybullet_data.getDataPath())
            self.world = pybullet.loadURDF("plane.urdf")
        pybullet.changeDynamics(self.world, -1, lateralFriction=1.0)
        pybullet.configureDebugVisualizer(pybullet.COV_ENABLE_GUI, 0)
        self.robot = pybullet.loadURDF(self.urdf_path,  useFixedBase=False, useMaximalCoordinates=False)
        self.init_state(q0)
        self.t = t0
        self.q = q0
        num_joints = len(self.get_joint_id_map())
        self.v = np.zeros(6+num_joints)
        self.u = np.zeros(num_joints)
        if self.camera_settings is not None:
            pybullet.resetDebugVisualizerCamera(self.camera_settings.camera_distance,
                                                self.camera_settings.camera_yaw,
                                                self.camera_settings.camera_pitch,
                                                self.camera_settings.camera_target_pos)

    def step_simulation(self, u: np.ndarray):
        base_lin_vel_world, _ = self.get_body_world_velocity()
        self.base_lin_vel_world_prev = base_lin_vel_world
        self.u = u
        self.apply_control_input(self.u)
        pybullet.stepSimulation()
        time.sleep(self.time_step)
        self.t += self.time_step

    def disconnect(self):
        pybullet.disconnect()

    def print_joint_info(self):
        pybullet.connect(pybullet.DIRECT)
        robot = pybullet.loadURDF(self.urdf_path, 
                                  useFixedBase=False, 
                                  useMaximalCoordinates=False)
        nJoints = pybullet.getNumJoints(robot)
        for j in range(nJoints):
            info = pybullet.getJointInfo(robot, j)
            print(info)
        pybullet.disconnect()