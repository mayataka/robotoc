#include "idocp/utils/trajectory_viewer.hpp"

#include <thread>
#include <chrono>
#include <cmath>

#include <boost/filesystem.hpp>

#include "pinocchio/gepetto/viewer.hpp"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/utils.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "gepetto/viewer/corba/client.hh"
#include "gepetto/viewer/corba/conversions.hh"


namespace idocp {

TrajectoryViewer::TrajectoryViewer(const std::string& path_to_urdf,
                                   const std::string& path_to_pkg)
  : force_radius_(0.015),
    force_length_(0.5),
    force_scale_(0.75),
    friction_cone_scale_(0.15) {
  pinocchio::urdf::buildModel(boost::filesystem::absolute(path_to_urdf).string(), 
                              model_);
  pinocchio::urdf::buildGeom(model_, 
                             boost::filesystem::absolute(path_to_urdf).string(), 
                             pinocchio::VISUAL, vmodel_, 
                             boost::filesystem::absolute(path_to_pkg).string());
  setForceProperties();
  setFrictionConeProperties();
  setCameraTransformDefault();
}


TrajectoryViewer::TrajectoryViewer(const std::string& path_to_urdf)
  : force_radius_(0.015),
    force_length_(0.5),
    force_scale_(0.75),
    friction_cone_scale_(0.15) {
  std::string path_to_pkg = boost::filesystem::absolute(path_to_urdf).string();
  for (int i=0; i<3; ++i) {
    const int file_name_length 
        = boost::filesystem::absolute(path_to_pkg).filename().string().size();
    for (int j=0; j<file_name_length; ++j) path_to_pkg.pop_back();
    path_to_pkg.pop_back();
  }
  pinocchio::urdf::buildModel(boost::filesystem::absolute(path_to_urdf).string(), 
                              model_);
  pinocchio::urdf::buildGeom(model_, 
                             boost::filesystem::absolute(path_to_urdf).string(), 
                             pinocchio::VISUAL, vmodel_, path_to_pkg);
  setForceProperties();
  setFrictionConeProperties();
  setCameraTransformDefault();
}


void TrajectoryViewer::setForceProperties() {
  force_radius_ = 0.015;
  force_length_ = 0.5;
  force_scale_ = 0.75;
  force_color_[0] = 1.0;
  force_color_[1] = 1.0;
  force_color_[2] = 0.0;
  force_color_[3] = 1.0;
}


void TrajectoryViewer::setFrictionConeProperties() {
  friction_cone_scale_ = 0.15;
  cone_color_[0] = 0.3;
  cone_color_[1] = 0.3;
  cone_color_[2] = 0.7;
  cone_color_[3] = 0.7;
  x_axis_ << 1, 0, 0;
}


void TrajectoryViewer::display(const std::vector<Eigen::VectorXd>& q_traj, 
                               const double sampling_period_in_sec) {
  pinocchio::gepetto::Viewer viewer(model_, &vmodel_, NULL);
  const bool success = viewer.initViewer("idocp::TrajectoryViewer");
  if(!success) {
    std::cout << "Failed in connecting to CORBA!!" << std::endl;
  }

  auto gui = gepetto::viewer::corba::gui();
  const auto window_id = gui->createWindow("idocp::TrajectoryViewer");
  gui->createScene("hpp-gui");
  gui->addSceneToWindow("hpp-gui", window_id);
  viewer.loadViewerModel("hpp-gui");
  for (int i=0; i<model_.nframes; ++i) {
    viewer.addFrame(i);
  }

  if (!gui->nodeExists("hpp-gui/floor")) {
    gui->addFloor("hpp-gui/floor");
  }
  gui->setColor("hpp-gui/floor", gepetto::viewer::corba::grey);
  gui->setBackgroundColor1("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setBackgroundColor2("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setLightingMode("hpp-gui/floor", "OFF");

  setCameraTransform();

  for (const auto& q : q_traj) {
    viewer.display(q);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}


void TrajectoryViewer::display(Robot& robot, 
                               const std::vector<Eigen::VectorXd>& q_traj, 
                               const std::vector<Eigen::VectorXd>& f_traj, 
                               const double sampling_period_in_sec,
                               const double mu) {
  pinocchio::gepetto::Viewer viewer(model_, &vmodel_, NULL);
  const bool success = viewer.initViewer("idocp::TrajectoryViewer");
  if(!success) {
    std::cout << "Failed in connecting to CORBA!!" << std::endl;
  }

  auto gui = gepetto::viewer::corba::gui();
  const auto window_id = gui->createWindow("idocp::TrajectoryViewer");
  gui->createScene("hpp-gui");
  gui->addSceneToWindow("hpp-gui", window_id);
  viewer.loadViewerModel("hpp-gui");
  for (int i=0; i<model_.nframes; ++i) {
    viewer.addFrame(i);
  }

  if (!gui->nodeExists("hpp-gui/floor")) {
    gui->addFloor("hpp-gui/floor");
  }
  gui->setColor("hpp-gui/floor", gepetto::viewer::corba::grey);
  gui->setBackgroundColor1("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setBackgroundColor2("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setLightingMode("hpp-gui/floor", "OFF");

  if (mu > 0) {
    gui->createGroup("hpp-gui/friction_cones");
    Eigen::Vector3d cone;
    cone << mu, mu, 1.0;
    for (int j=0; j<robot.maxPointContacts(); ++j) {
        gui->createGroup(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)).c_str());
        gui->addCurve(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/vertex").c_str(), 
                      gepetto::viewer::corba::positionSeq({ {0., 0., 0., }, 
                                                            { cone[0],  cone[1],  cone[2]}, 
                                                            {-cone[0],  cone[1],  cone[2]}, 
                                                            {-cone[0], -cone[1],  cone[2]}, 
                                                            { cone[0], -cone[1],  cone[2]}, 
                                                            { cone[0],  cone[1],  cone[2]}, }), cone_color_);
        gui->setCurveMode(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/vertex").c_str(), "TRIANGLE_FAN");
        gepetto::corbaserver::Position p0, p;
        p0[0] = 0.; p0[1] = 0.; p0[2] = 0.;
        p[0] =  cone[0]; p[1] =  cone[1]; p[2] = cone[2];
        gui->addLine(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line1").c_str(), p0, p, cone_color_);
        p[0] = -cone[0]; p[1] =  cone[1]; p[2] = cone[2];
        gui->addLine(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line2").c_str(), p0, p, cone_color_);
        p[0] = -cone[0]; p[1] = -cone[1]; p[2] = cone[2];
        gui->addLine(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line3").c_str(), p0, p, cone_color_);
        p[0] =  cone[0]; p[1] = -cone[1]; p[2] = cone[2];
        gui->addLine(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line4").c_str(), p0, p, cone_color_);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)).c_str(), "Alpha", 0);
        gepetto::corbaserver::Position cone_scale;
        cone_scale[0] = friction_cone_scale_; cone_scale[1] = friction_cone_scale_; cone_scale[2] = friction_cone_scale_;
        gui->setScale(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/vertex").c_str(), cone_scale);
        gui->setScale(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line1").c_str(), cone_scale);
        gui->setScale(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line2").c_str(), cone_scale);
        gui->setScale(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line3").c_str(), cone_scale);
        gui->setScale(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line4").c_str(), cone_scale);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/vertex").c_str(), "Alpha", 0);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line1").c_str(), "Alpha", 0.2);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line2").c_str(), "Alpha", 0.2);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line3").c_str(), "Alpha", 0.2);
        gui->setFloatProperty(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)+"/line4").c_str(), "Alpha", 0.2);
    }
  }

  gui->createGroup("hpp-gui/contact_forces");
  for (int j=0; j<robot.maxPointContacts(); ++j) {
    gui->addArrow(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), force_radius_, force_length_, force_color_);
    gui->setFloatProperty(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), "Alpha", 1.0);
    gui->setVisibility(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), "ALWAYS_ON_TOP");
  }

  setCameraTransform();

  for (int i=0; i<f_traj.size(); ++i) {
    robot.updateFrameKinematics(q_traj[i]);
    for (int j=0; j<robot.maxPointContacts(); ++j) {
      const Eigen::Vector3d& f = f_traj[i].template segment<3>(3*j);
      gepetto::corbaserver::Position f_scale;
      f_scale[0] = force_scale_ * std::sqrt(f.norm() / robot.totalWeight());
      f_scale[1] = 1.;
      f_scale[2] = 1.;
      gui->setVector3Property(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), "Scale", f_scale);
      const Eigen::Matrix3d R = rotationMatrix(f);
      const pinocchio::SE3::Quaternion quat = pinocchio::SE3::Quaternion(R);
      gepetto::corbaserver::Transform pose;
      const int contact_frame = robot.contactFrames()[j];
      pose[0] = robot.framePosition(contact_frame)[0];
      pose[1] = robot.framePosition(contact_frame)[1];
      pose[2] = robot.framePosition(contact_frame)[2];
      pose[3] = quat.x();
      pose[4] = quat.y();
      pose[5] = quat.z();
      pose[6] = quat.w();
      gui->applyConfiguration(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), pose);
      if (mu > 0) {
        pose[0] = robot.framePosition(contact_frame)[0];
        pose[1] = robot.framePosition(contact_frame)[1];
        pose[2] = robot.framePosition(contact_frame)[2];
        pose[3] = 0.;
        pose[4] = 0.;
        pose[5] = 0.;
        pose[6] = 1.;
        gui->applyConfiguration(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)).c_str(), pose);
      }
      if (f.norm() > 0) {
        if (mu > 0) {
          gui->setVisibility(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)).c_str(), "ON");
        }
        gui->setVisibility(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), "ON");
      }
      else {
        if (mu > 0) {
          gui->setVisibility(("hpp-gui/friction_cones/friction_cone_"+std::to_string(j)).c_str(), "OFF");
        }
        gui->setVisibility(("hpp-gui/contact_forces/contact_force_"+std::to_string(j)).c_str(), "OFF");
      }
      gui->refresh();
    }
    viewer.display(q_traj[i]);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}


void TrajectoryViewer::setCameraTransform(const Eigen::Vector3d& camera_pos,
                                          const Eigen::Vector4d& camera_quat) {
  camera_pos_ = camera_pos;
  camera_quat_ = camera_quat;
}


void TrajectoryViewer::setCameraTransformDefault() {
  camera_pos_ << 2.2, -3.5, 1.13;
  camera_quat_ << 0.60612, 0.166663, 0.19261, 0.753487;
}


void TrajectoryViewer::printCurrentCameraTransform() const {
  auto gui = gepetto::viewer::corba::gui();
  const auto window_id = gui->getWindowID("idocp::TrajectoryViewer");
  const auto camera = gui->getCameraTransform(window_id);
  std::cout << "Current camera transform is: [";
  for (int i=0; i<6; ++i) {
    std::cout << camera[i] << ", ";
  }
  std::cout << camera[6] << "]" << std::endl;
}


void TrajectoryViewer::setCameraTransform() {
  auto gui = gepetto::viewer::corba::gui();
  const auto window_id = gui->getWindowID("idocp::TrajectoryViewer");
  std::vector<float> pose(7);
  pose[0] = camera_pos_.coeff(0);
  pose[1] = camera_pos_.coeff(1);
  pose[2] = camera_pos_.coeff(2);
  pose[3] = camera_quat_.normalized().coeff(0);
  pose[4] = camera_quat_.normalized().coeff(1);
  pose[5] = camera_quat_.normalized().coeff(2);
  pose[6] = camera_quat_.normalized().coeff(3);
  gui->setCameraTransform(window_id, pose.data());
}

} // namespace idocp