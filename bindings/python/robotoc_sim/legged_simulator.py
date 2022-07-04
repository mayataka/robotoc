import pybullet
import math
import time
import abc
import numpy as np
import os


class LeggedSimulator(metaclass=abc.ABCMeta):
    def __init__(self, path_to_urdf, time_step, start_time, end_time):
        self.path_to_urdf = path_to_urdf
        self.time_step = time_step
        self.start_time = start_time
        self.end_time = end_time
        self.calib_camera = False
        self.camera_distance = 0.0
        self.camera_yaw = 0.0
        self.camera_pitch = 0.0
        self.camera_target_pos = [0., 0., 0.]
        self.print_items = []
        self.terrain_urdf = os.path.join(os.path.dirname(__file__), "rsc/terrain.urdf")
        self.height_range = None
        self.num_heightfield_rows = None
        self.num_heightfield_columns = None
        self.log_dir = None
        self.q_log = None
        self.v_log = None
        self.u_log = None
        self.t_log = None
        self.kkt_log = None
        self.sim_name = None

    def set_sim_settings(self, time_step, start_time, end_time):
        self.time_step  = time_step
        self.start_time = start_time
        self.end_time   = end_time

    def set_urdf(self, path_to_urdf):
        self.path_to_urdf = path_to_urdf

    def set_camera(self, camera_distance, camera_yaw, camera_pitch, camera_target_pos):
        self.calib_camera = True
        self.camera_distance = camera_distance
        self.camera_yaw = camera_yaw
        self.camera_pitch = camera_pitch
        self.camera_target_pos = camera_target_pos

    def get_body_local_velocity(pybullet_robot):
        basePos, baseOrn = pybullet.getBasePositionAndOrientation(pybullet_robot)
        R = np.reshape(pybullet.getMatrixFromQuaternion(baseOrn), [3, 3]) 
        baseVel, baseAngVel = pybullet.getBaseVelocity(pybullet_robot)
        baseVel = R.T @ np.array(baseVel)
        baseAngVel = R.T @ np.array(baseAngVel)
        return baseVel, baseAngVel

    @abc.abstractmethod
    def get_state_from_pybullet(self, pybullet_robot):
        return NotImplementedError()

    @abc.abstractmethod
    def apply_control_input_to_pybullet(self, pybullet_robot, u):
        return NotImplementedError()

    @abc.abstractmethod
    def init_pybullet_robot(self, pybullet_robot, q):
        return NotImplementedError()

    def print_joint_info(self):
        pybullet.connect(pybullet.DIRECT)
        robot = pybullet.loadURDF(self.path_to_urdf, 
                                  useFixedBase=False, 
                                  useMaximalCoordinates=False)
        nJoints = pybullet.getNumJoints(robot)
        for j in range(nJoints):
            info = pybullet.getJointInfo(robot, j)
            print(info)
        pybullet.disconnect()

    def add_print_item(self, item):
        self.print_items.append(item)

    def run_simulation(self, mpc, q0, v0, feedback_delay=False, terrain=False, 
                       verbose=False, log=False, record=False, 
                       sim_name='mpc_sim'):
        pybullet.connect(pybullet.GUI)
        pybullet.setGravity(0, 0, -9.81)
        pybullet.setTimeStep(self.time_step)
        if terrain:
            if self.height_range is not None \
                and self.num_heightfield_rows is not None \
                and self.num_heightfield_columns is not None:

                terrain_shape = self._create_random_terrain_shape(self.height_range, 
                                                                  self.num_heightfield_rows, 
                                                                  self.num_heightfield_columns)
                terrain = pybullet.createMultiBody(0, terrain_shape)
                pybullet.changeVisualShape(terrain, -1, rgbaColor=[1.0, 1.0, 1.0, 1.0])
                pybullet.resetBasePositionAndOrientation(terrain, [0., 0., 0.], [0., 0., 0., 1.])
                ray_test = pybullet.rayTest([q0[0], q0[1], self.height_range], 
                                            [q0[0], q0[1], -self.height_range])
                ray_hit_position = ray_test[0][3]
                pybullet.resetBasePositionAndOrientation(terrain, [0., 0., -ray_hit_position[2]], [0., 0., 0., 1.])
                pybullet.changeDynamics(terrain, -1, lateralFriction=1.0)
            else:
                terrain = pybullet.loadURDF(fileName=self.terrain_urdf)
        else:
            import pybullet_data
            pybullet.setAdditionalSearchPath(pybullet_data.getDataPath())
            plane = pybullet.loadURDF("plane.urdf")
        robot = pybullet.loadURDF(self.path_to_urdf,  
                                  useFixedBase=False, 
                                  useMaximalCoordinates=False)

        self.init_pybullet_robot(robot, q0)
        t = self.start_time
        q = q0
        v = v0
        u = mpc.get_initial_control_input().copy()
        sim_time = self.end_time - self.start_time
        sim_steps = math.floor(sim_time/self.time_step)

        if log:
            log_dir = os.path.join(os.getcwd(), sim_name+"_log")
            self.log_dir = log_dir
            os.makedirs(log_dir, exist_ok=True)
            q_log = open(os.path.join(log_dir, "q.log"), mode='w')
            v_log = open(os.path.join(log_dir, "v.log"), mode='w')
            u_log = open(os.path.join(log_dir, "u.log"), mode='w')
            t_log = open(os.path.join(log_dir, "t.log"), mode='w')
            kkt_log = open(os.path.join(log_dir, "kkt.log"), mode='w')

        if self.calib_camera:
            pybullet.resetDebugVisualizerCamera(self.camera_distance,
                                                self.camera_yaw,
                                                self.camera_pitch,
                                                self.camera_target_pos)
        pybullet.configureDebugVisualizer(pybullet.COV_ENABLE_GUI, 0)

        if record:
            pybullet.startStateLogging(pybullet.STATE_LOGGING_VIDEO_MP4, 
                                       sim_name+".mp4")

        for i in range(sim_steps):
            q, v = self.get_state_from_pybullet(robot)
            assert q.shape[0] == q0.shape[0]
            assert v.shape[0] == v0.shape[0]
            if verbose:
                print('t = {:.6g}:'.format(t))
            if feedback_delay:
                u = mpc.get_initial_control_input().copy()
            mpc.update_solution(t, self.time_step, q, v)
            kkt_error = mpc.KKT_error(t, q, v) 
            if verbose:
                print('KKT error = {:.6g}'.format(kkt_error))
                print('')
                if self.print_items:
                    for e in self.print_items:
                        print(e)
            if not feedback_delay:
                u = mpc.get_initial_control_input().copy()
            self.apply_control_input_to_pybullet(robot, u)
            pybullet.stepSimulation()
            time.sleep(self.time_step)
            if log:
                np.savetxt(q_log, [q])
                np.savetxt(v_log, [v])
                np.savetxt(u_log, [u])
                np.savetxt(t_log, np.array([t]))
                np.savetxt(kkt_log, np.array([kkt_error]))
            t = t + self.time_step

        if log:
            q_log.close()
            v_log.close()
            u_log.close()
            t_log.close()
            kkt_log.close()
            self.q_log = os.path.abspath(os.path.join(log_dir, "q.log"))
            self.v_log = os.path.abspath(os.path.join(log_dir, "v.log"))
            self.u_log = os.path.abspath(os.path.join(log_dir, "u.log"))
            self.t_log = os.path.abspath(os.path.join(log_dir, "t.log"))
            self.kkt_log = os.path.abspath(os.path.join(log_dir, "kkt.log"))
            self.log_dir = os.path.abspath(log_dir)
            print('Logs are saved at ' + self.log_dir)
        
        self.sim_name = sim_name


    def disconnect(self):
        pybullet.disconnect()


    def set_terrain_shape(self, height_range=1.0, 
                          num_heightfield_rows=100,
                          num_heightfield_columns=100):
        self.height_range = height_range
        self.num_heightfield_rows = num_heightfield_rows
        self.num_heightfield_columns = num_heightfield_columns


    def _create_random_terrain_shape(self, height_range, 
                                     num_heightfield_rows, 
                                     num_heightfield_columns):
        heightfield_data = np.zeros(shape=[num_heightfield_columns, num_heightfield_rows], dtype=np.float)
        for i in range(int(num_heightfield_columns/2)):
            for j in range(int(num_heightfield_rows)):
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
            heightfield_data = (2.0 * height_range / (max-min)) * heightfield_data 
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