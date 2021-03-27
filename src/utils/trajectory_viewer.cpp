#include "idocp/utils/trajectory_viewer.hpp"

#include <thread>
#include <chrono>
#include <cmath>

#include <boost/filesystem.hpp>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/utils.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/gepetto/viewer.hpp"
#include "gepetto/viewer/corba/conversions.hh"


namespace idocp {

TrajectoryViewer::TrajectoryViewer(
    const std::string& description_pkg_serach_path, 
    const std::string& path_to_urdf)
  : description_pkg_serach_path_(description_pkg_serach_path),
    path_to_urdf_(path_to_urdf) {
  cone_color_[0] = 0.3;
  cone_color_[1] = 0.3;
  cone_color_[2] = 0.7;
  cone_color_[3] = 0.7;
  force_radius_ = 0.015;
  force_length_ = 0.5;
  friction_cone_scale_ = 0.15;
  x_axis_ << 1, 0, 0;
}


TrajectoryViewer::TrajectoryViewer(const std::string& path_to_urdf)
  : description_pkg_serach_path_(boost::filesystem::absolute("../"+path_to_urdf).string()),
    path_to_urdf_(path_to_urdf) {
  cone_color_[0] = 0.3;
  cone_color_[1] = 0.3;
  cone_color_[2] = 0.7;
  cone_color_[3] = 0.7;
  force_radius_ = 0.015;
  force_length_ = 0.5;
  friction_cone_scale_ = 0.15;
  x_axis_ << 1, 0, 0;
}


void TrajectoryViewer::display(const std::vector<Eigen::VectorXd>& q_traj, 
                               const double sampling_period_in_sec) {
  pinocchio::Model model;
  pinocchio::urdf::buildModel(path_to_urdf_, model);
  pinocchio::GeometryModel vmodel;
  pinocchio::urdf::buildGeom(model, path_to_urdf_, pinocchio::VISUAL, vmodel, 
                             description_pkg_serach_path_);

  pinocchio::gepetto::Viewer viewer(model, &vmodel, NULL);
  const bool success = viewer.initViewer("idocp::TrajectoryViewer");
  if(!success) {
    std::cout << "Failed in connecting to CORBA!!" << std::endl;
  }
  viewer.loadViewerModel();

  for (int i=0; i<model.nframes; ++i) {
    viewer.addFrame(i);
  }
  // viewer.addFrame(model.nframes-1);

  auto gui = gepetto::viewer::corba::gui();
  gui->setColor("hpp-gui/floor", gepetto::viewer::corba::grey);
  gui->setBackgroundColor1("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setBackgroundColor2("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setLightingMode("hpp-gui/floor", "OFF");

  for (const auto& q : q_traj) {
    viewer.display(q);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}


void TrajectoryViewer::display(Robot& robot, const double mu,
                               const std::vector<Eigen::VectorXd>& q_traj, 
                               const std::vector<Eigen::VectorXd>& f_traj, 
                               const double sampling_period_in_sec) {
  pinocchio::Model model;
  pinocchio::urdf::buildModel(path_to_urdf_, model);
  pinocchio::GeometryModel vmodel;
  pinocchio::urdf::buildGeom(model, path_to_urdf_, pinocchio::VISUAL, vmodel, 
                             description_pkg_serach_path_);

  pinocchio::gepetto::Viewer viewer(model, &vmodel, NULL);
  const bool success = viewer.initViewer("idocp::TrajectoryViewer");
  if(!success) {
    std::cout << "Failed in connecting to CORBA!!" << std::endl;
  }
  viewer.loadViewerModel();

  for (int i=0; i<model.nframes; ++i) {
    viewer.addFrame(i);
  }
  // viewer.addFrame(model.nframes-1);

  auto gui = gepetto::viewer::corba::gui();
  gui->setColor("hpp-gui/floor", gepetto::viewer::corba::grey);
  gui->setBackgroundColor1("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setBackgroundColor2("idocp::TrajectoryViewer", gepetto::viewer::corba::white);
  gui->setLightingMode("hpp-gui/floor", "OFF");
  Eigen::Vector3d cone;
  cone << mu, mu, 1.0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    gui->createGroup(("hpp-gui/contact_force"+std::to_string(i)).c_str());
    gui->createGroup(("hpp-gui/friction_cone"+std::to_string(i)).c_str());
    gui->addCurve(("hpp-gui/friction_cone"+std::to_string(i)+"/vertex").c_str(), 
                  gepetto::viewer::corba::positionSeq({ {0., 0., 0., }, 
                                                        { cone[0],  cone[1],  cone[2]}, 
                                                        {-cone[0],  cone[1],  cone[2]}, 
                                                        {-cone[0], -cone[1],  cone[2]}, 
                                                        { cone[0], -cone[1],  cone[2]}, 
                                                        { cone[0],  cone[1],  cone[2]}, }), cone_color_);
    gui->setCurveMode(("hpp-gui/friction_cone"+std::to_string(i)+"/vertex").c_str(), "TRIANGLE_FAN");
    gepetto::corbaserver::Position p0, p;
    p0[0] = 0.; p0[1] = 0.; p0[2] = 0.;
    p[0] =  cone[0]; p[1] =  cone[1]; p[2] = cone[2];
    gui->addLine(("hpp-gui/friction_cone"+std::to_string(i)+"/line1").c_str(), p0, p, cone_color_);
    p[0] = -cone[0]; p[1] =  cone[1]; p[2] = cone[2];
    gui->addLine(("hpp-gui/friction_cone"+std::to_string(i)+"/line2").c_str(), p0, p, cone_color_);
    p[0] = -cone[0]; p[1] = -cone[1]; p[2] = cone[2];
    gui->addLine(("hpp-gui/friction_cone"+std::to_string(i)+"/line3").c_str(), p0, p, cone_color_);
    p[0] =  cone[0]; p[1] = -cone[1]; p[2] = cone[2];
    gui->addLine(("hpp-gui/friction_cone"+std::to_string(i)+"/line4").c_str(), p0, p, cone_color_);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)).c_str(), "Alpha", 0);
    gepetto::corbaserver::Position cone_scale;
    cone_scale[0] = friction_cone_scale_; cone_scale[1] = friction_cone_scale_; cone_scale[2] = friction_cone_scale_;
    gui->setScale(("hpp-gui/friction_cone"+std::to_string(i)+"/vertex").c_str(), cone_scale);
    gui->setScale(("hpp-gui/friction_cone"+std::to_string(i)+"/line1").c_str(), cone_scale);
    gui->setScale(("hpp-gui/friction_cone"+std::to_string(i)+"/line2").c_str(), cone_scale);
    gui->setScale(("hpp-gui/friction_cone"+std::to_string(i)+"/line3").c_str(), cone_scale);
    gui->setScale(("hpp-gui/friction_cone"+std::to_string(i)+"/line4").c_str(), cone_scale);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)+"/vertex").c_str(), "Alpha", 0);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)+"/line1").c_str(), "Alpha", 0.2);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)+"/line2").c_str(), "Alpha", 0.2);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)+"/line3").c_str(), "Alpha", 0.2);
    gui->setFloatProperty(("hpp-gui/friction_cone"+std::to_string(i)+"/line4").c_str(), "Alpha", 0.2);
    gui->addArrow(("hpp-gui/contact_force"+std::to_string(i)).c_str(), force_radius_, force_length_, gepetto::viewer::corba::yellow);
    gui->setFloatProperty(("hpp-gui/contact_force"+std::to_string(i)).c_str(), "Alpha", 1.0);
    gui->setVisibility(("hpp-gui/contact_force"+std::to_string(i)).c_str(), "ALWAYS_ON_TOP");
  }

  for (int i=0; i<q_traj.size(); ++i) {
    robot.updateFrameKinematics(q_traj[i]);
    for (int j=0; j<robot.maxPointContacts(); ++j) {
      const int contact_frame = robot.contactFrames()[j];
      const Eigen::Vector3d f = robot.frameRotation(contact_frame) * f_traj[i].template segment<3>(3*j);
      gepetto::corbaserver::Position f_scale;
      f_scale[0] = std::sqrt(f.norm() / robot.totalWeight());
      f_scale[1] = 1.;
      f_scale[2] = 1.;
      gui->setVector3Property(("hpp-gui/contact_force"+std::to_string(j)).c_str(), "Scale", f_scale);
      const Eigen::Matrix3d R = rotationMatrix(x_axis_, f);
      const pinocchio::SE3::Quaternion quat = pinocchio::SE3::Quaternion(R);
      gepetto::corbaserver::Transform pose;
      pose[0] = f[0];
      pose[1] = f[1];
      pose[2] = f[2];
      pose[3] = quat.x();
      pose[4] = quat.y();
      pose[5] = quat.z();
      pose[6] = quat.w();
      gui->applyConfiguration(("hpp-gui/contact_force"+std::to_string(j)).c_str(), pose);
      pose[0] = robot.framePosition(contact_frame)[0];
      pose[1] = robot.framePosition(contact_frame)[1];
      pose[2] = robot.framePosition(contact_frame)[2];
      pose[3] = 0.;
      pose[4] = 0.;
      pose[5] = 0.;
      pose[6] = 1.;
      gui->applyConfiguration(("hpp-gui/friction_cone"+std::to_string(j)).c_str(), pose);
      if (f.norm() > 0) {
        gui->setVisibility(("hpp-gui/friction_cone"+std::to_string(j)).c_str(), "ON");
        gui->setVisibility(("hpp-gui/contact_force"+std::to_string(j)).c_str(), "ON");
      }
      else {
        gui->setVisibility(("hpp-gui/friction_cone"+std::to_string(j)).c_str(), "OFF");
        gui->setVisibility(("hpp-gui/contact_force"+std::to_string(j)).c_str(), "OFF");
      }
      gui->refresh();
    }
    viewer.display(q_traj[i]);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}

} // namespace idocp