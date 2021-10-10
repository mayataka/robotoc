import idocp
import pinocchio
from pinocchio.robot_wrapper import RobotWrapper
from pinocchio.visualize import GepettoVisualizer, MeshcatVisualizer
from os.path import abspath, dirname, join
import numpy as np
import math
import time
import meshcat.transformations 

class TrajectoryViewer:
    def __init__(self, path_to_urdf, path_to_pkg=None, 
                 base_joint_type=idocp.BaseJointType.FixedBase, 
                 viewer_type='gepetto'):
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

        self.viewer_type = viewer_type
        if viewer_type == 'gepetto':
            self.viewer = GepettoVisualizer(self.robot.model, 
                                            self.robot.collision_model,
                                            self.robot.visual_model)
            self.window_name = 'idocp.TrajectoryViewer'
            self.camera_pos = [2.2, -3.5, 1.13] 
            self.camera_angle = [0.60612, 0.166663, 0.19261, 0.753487] 
        elif viewer_type == 'meshcat':
            self.viewer = MeshcatVisualizer(self.robot.model, 
                                            self.robot.collision_model, 
                                            self.robot.visual_model)
            self.camera_tf = meshcat.transformations.translation_matrix([0.8, -2.5, -0.2]) 
            self.zoom = 3.0
        else:
            print('Please choose viewer_type from "gepetto" or "meshcat"!')
            return NotImplementedError()


    def set_contact_info(self, contact_frames, mu):
        self.display_contact = True
        self.contact_frames = contact_frames
        self.mu = mu


    def set_camera_transform_gepetto(self, camera_pos=None, camera_angle=None):
        if camera_pos is not None:
            self.camera_pos = camera_pos
        if camera_angle is not None:
            self.camera_angle = camera_angle

    def set_camera_transform_meshcat(self, camera_tf_vec=None, zoom=None):
        if camera_tf_vec is not None:
            self.camera_tf = meshcat.transformations.translation_matrix(camera_tf_vec) 
        if zoom is not None:
            self.zoom = zoom


    def display(self, dt, q_traj, f_traj=None):
        if self.viewer_type == 'gepetto':
            self.display_gepetto(dt, q_traj, f_traj)
        elif self.viewer_type == 'meshcat':
            self.display_meshcat(dt, q_traj)


    def display_gepetto(self, dt, q_traj, f_traj):
        self.robot.setVisualizer(self.viewer)
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
                    f_scale = [math.sqrt(self.force_scale*np.linalg.norm(fi)/robot.total_weight()), 1.0, 1.0]
                    gui.setVector3Property('world/contact_forces/contact_force_'+str(i), "Scale", f_scale)
                    fpos = robot.frame_position(self.contact_frames[i])
                    quat = pinocchio.Quaternion(self.x_axis, fi).normalized()
                    pose = np.concatenate((fpos, np.array([quat.x, quat.y, quat.z, quat.w])), axis=None).tolist()
                    gui.applyConfiguration('world/contact_forces/contact_force_'+str(i), pose)
                    if self.mu > 0.0:
                        pose = np.concatenate((fpos, np.array([0., 0., 0., 1.])), axis=None).tolist()
                        gui.applyConfiguration('world/friction_cones/friction_cone_'+str(i), pose)
                    if np.linalg.norm(fi) > 0.0 and self.mu > 0.0:
                        gui.setVisibility('world/friction_cones/friction_cone_'+str(i), 'ON')
                        gui.setVisibility('world/contact_forces/contact_force_'+str(i), 'ON')
                    else:
                        gui.setVisibility('world/friction_cones/friction_cone_'+str(i), 'OFF')
                        gui.setVisibility('world/contact_forces/contact_force_'+str(i), 'OFF')
                    gui.refresh()
                self.robot.display(q)
                time.sleep(sleep_time)
        else:
            for q in q_traj:
                self.robot.display(q)
                time.sleep(sleep_time)


    def display_meshcat(self, dt, q_traj, open=True):
        self.robot.setVisualizer(self.viewer)
        self.robot.initViewer(open=open)
        self.robot.loadViewerModel(rootNodeName='idocp.TrajectoryViewer')
        self.viewer.viewer["/Cameras/default"].set_transform(self.camera_tf)
        self.viewer.viewer["/Cameras/default/rotated/<object>"].set_property("zoom", self.zoom)
        # self.viewer.viewer["/Cameras/default"].set_transform(tf.translation_matrix([0.3, -2.5, -0.2]))
        # self.viewer.viewer["/Cameras/default/rotated/<object>"].set_property("zoom", 5.0)

        sleep_time = dt / self.play_speed
        for q in q_traj:
            self.robot.display(q)
            time.sleep(sleep_time)