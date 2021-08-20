import idocp
import pinocchio
from pinocchio.robot_wrapper import RobotWrapper
from pinocchio.utils import skew
from pinocchio.visualize import GepettoVisualizer, MeshcatVisualizer
from os.path import abspath, dirname, join
import numpy as np
import time


class TrajectoryViewer:
    def __init__(self, path_to_urdf, path_to_pkg=None, 
                 base_joint_type=idocp.BaseJointType.FixedBase):
        self.path_to_urdf = abspath(path_to_urdf)
        if path_to_pkg is None:
            path_to_pkg = join(dirname(self.path_to_urdf), '../..')
        if base_joint_type == idocp.BaseJointType.FloatingBase:
            self.robot= RobotWrapper.BuildFromURDF(self.path_to_urdf, path_to_pkg,
                                                   pinocchio.JointModelFreeFlyer())
        else:
            self.robot= RobotWrapper.BuildFromURDF(self.path_to_urdf, path_to_pkg)
        self.path_to_urdf = path_to_urdf
        self.base_joint_type = base_joint_type

        self.play_speed = 1.0

        self.window_name = 'idocp.TrajectoryViewer'
        self.camera_pos=[2.2, -3.5, 1.13]
        self.camera_angle=[0.60612, 0.166663, 0.19261, 0.753487]

        self.force_radius = 0.015
        self.force_length = 0.5
        self.force_scale = 0.75
        self.force_color = [1.0, 1.0, 0.0, 1.0]

        self.display_contact = False
        self.contact_frames = []
        self.mu = 0.
        self.cone_scale = [0.15, 0.15, 0.15]
        self.cone_color = [0.3, 0.3, 0.7, 0.7]
        self.x_axis = np.array([1.0, 0.0, 0.0])


    def set_contact_info(self, contact_frames, mu):
        self.display_contact = True
        self.contact_frames = contact_frames
        self.mu = mu


    def set_camera_transform(self, camera_pos, camera_angle):
        self.camera_pos = camera_pos
        self.camera_angle = camera_angle


    def display(self, dt, q_traj, f_traj=None, viewer='gepetto'):
        if viewer == 'gepetto':
            self.display_gepetto(dt, q_traj, f_traj)
        elif viewer == 'meshcat':
            self.display_meshcat(dt, q_traj)


    def rotation_matrix(self, vec):
            nrm = np.linalg.norm(vec)
            vec_nrmrzd = np.array(vec/nrm) if nrm > 0.0 else np.zeros(3)
            axis_cross_vec = np.cross(self.x_axis, vec_nrmrzd)
            cross_norm = np.linalg.norm(axis_cross_vec)
            if cross_norm == 0.:
                return np.eye(3, 3)
            else:
                inner = np.dot(self.x_axis, vec_nrmrzd)
                skew_mat = pinocchio.skew(axis_cross_vec)
                return np.eye(3, 3) + skew_mat + skew_mat*skew_mat*((1.0-inner)/(cross_norm*cross_norm))


    def display_gepetto(self, dt, q_traj, f_traj):
        viz = GepettoVisualizer(self.robot.model, self.robot.collision_model,
                                self.robot.visual_model)
        self.robot.setVisualizer(viz)
        self.robot.initViewer(windowName=self.window_name, loadModel=False)
        self.robot.loadViewerModel(rootNodeName=self.window_name)
        gui = self.robot.viewer.gui
        # init the window
        window_id = gui.getWindowID(self.window_name)
        gui.setBackgroundColor1(window_id, [1., 1., 1., 1.])
        gui.setBackgroundColor2(window_id, [1., 1., 1., 1.])
        # init the floor
        floor_name = 'world/floor'
        gui.addFloor(floor_name)
        gui.setColor(floor_name, [0.7, 0.7, 0.7, 1.0])
        gui.setLightingMode(floor_name, 'OFF')
        # init contact forces and friction cones
        if f_traj is not None and self.display_contact:
            # create cones
            self.robot.viewer.gui.createGroup('world/friction_cones')
            cone = [self.mu, self.mu, 1.0]
            for i in range(len(self.contact_frames)):
                gui.createGroup('world/friction_cones/friction_cone_'+str(i))
                gui.addCurve('world/friction_cones/friction_cone_'+str(i)+'/vertex', 
                             [[0., 0., 0.], 
                             [ cone[0],  cone[1],  cone[2]],
                             [-cone[0],  cone[1],  cone[2]],
                             [-cone[0], -cone[1],  cone[2]],
                             [ cone[0], -cone[1],  cone[2]],
                             [ cone[0],  cone[1],  cone[2]]], self.cone_color)
                gui.setCurveMode('world/friction_cones/friction_cone_'+str(i)+'/vertex', 'TRIANGLE_FAN')
                gui.addLine('world/friction_cones/friction_cone_'+str(i)+'/line1', 
                            [0., 0., 0.], [ cone[0],  cone[1],  cone[2]], self.cone_color)
                gui.addLine('world/friction_cones/friction_cone_'+str(i)+'/line2', 
                            [0., 0., 0.], [-cone[0],  cone[1],  cone[2]], self.cone_color)
                gui.addLine('world/friction_cones/friction_cone_'+str(i)+'/line3', 
                            [0., 0., 0.], [-cone[0], -cone[1],  cone[2]], self.cone_color)
                gui.addLine('world/friction_cones/friction_cone_'+str(i)+'/line4', 
                            [0., 0., 0.], [ cone[0], -cone[1],  cone[2]], self.cone_color)
                gui.setScale('world/friction_cones/friction_cone_'+str(i)+'/vertex', self.cone_scale)
                gui.setScale('world/friction_cones/friction_cone_'+str(i)+'/line1', self.cone_scale)
                gui.setScale('world/friction_cones/friction_cone_'+str(i)+'/line2', self.cone_scale)
                gui.setScale('world/friction_cones/friction_cone_'+str(i)+'/line3', self.cone_scale)
                gui.setScale('world/friction_cones/friction_cone_'+str(i)+'/line4', self.cone_scale)
                gui.setFloatProperty('world/friction_cones/friction_cone_'+str(i)+'/vertex', 'Alpha', 0.)
                gui.setFloatProperty('world/friction_cones/friction_cone_'+str(i)+'/line1', 'Alpha', 0.2)
                gui.setFloatProperty('world/friction_cones/friction_cone_'+str(i)+'/line2', 'Alpha', 0.2)
                gui.setFloatProperty('world/friction_cones/friction_cone_'+str(i)+'/line3', 'Alpha', 0.2)
                gui.setFloatProperty('world/friction_cones/friction_cone_'+str(i)+'/line4', 'Alpha', 0.2)
            # create forces
            gui.createGroup('world/contact_forces')
            for i in range(len(self.contact_frames)):
                gui.addArrow('world/contact_forces/contact_force_'+str(i),
                             self.force_radius, self.force_length, self.force_color)
                gui.setFloatProperty('world/contact_forces/contact_force_'+str(i), 'Alpha', 1.0)
                gui.setVisibility('world/contact_forces/contact_force_'+str(i), 'ALWAYS_ON_TOP')
        # set camera and play speed
        sleep_time = dt / self.play_speed
        camera = self.camera_pos
        camera.extend(self.camera_angle)
        gui.setCameraTransform(self.robot.viz.windowID, camera)
        # display
        if f_traj is not None:
            robot = idocp.Robot(self.path_to_urdf, self.base_joint_type, 
                                self.contact_frames)
            for q, f in zip(q_traj, f_traj):
                robot.forward_kinematics(q)
                for i in range(len(self.contact_frames)):
                    fi = f[3*i:3*(i+1)]
                    f_scale = [(self.force_scale*np.linalg.norm(fi)/robot.total_weight()), 1.0, 1.0]
                    gui.setVector3Property('world/contact_forces/contact_force_'+str(i), "Scale", f_scale)
                    fpos = robot.frame_position(self.contact_frames[i])
                    quat = pinocchio.Quaternion(self.x_axis, fi).normalized()
                    pose = np.concatenate((fpos, np.array([quat.x, quat.y, quat.z, quat.w])), axis=None)
                    gui.applyConfiguration('world/contact_forces/contact_force_'+str(i), pose.tolist())
                    if self.mu > 0.0:
                        pose = np.concatenate((fpos, np.array([0., 0., 0., 1.])), axis=None)
                        gui.applyConfiguration('world/friction_cones/friction_cone_'+str(i), pose.tolist())
                    if np.linalg.norm(fi) > 0.0 and self.mu > 0.0:
                        gui.setVisibility('world/friction_cones/friction_cone_'+str(i), 'ON')
                        gui.setVisibility('world/contact_forces/contact_force_'+str(i), 'ON')
                    else:
                        gui.setVisibility('world/friction_cones/friction_cone_'+str(i), 'OFF')
                        gui.setVisibility('world/contact_forces/contact_force_'+str(i), 'OFF')
                    gui.refresh()
                viz.display(q)
                time.sleep(sleep_time)
        else:
            for q in q_traj:
                viz.display(q)
                time.sleep(sleep_time)


    def display_meshcat(self, dt, q_traj, open=True):
        viz = MeshcatVisualizer(self.robot.model, self.robot.collision_model, 
                                self.robot.visual_model)
        self.robot.setVisualizer(viz)
        self.robot.initViewer(open=open)
        self.robot.loadViewerModel(rootNodeName='idocp.TrajectoryViewer')
        sleep_time = dt / self.play_speed
        for q in q_traj:
            viz.display(q)
            time.sleep(sleep_time)