#include "idocp/utils/trajectory_viewer.hpp"

#include <thread>
#include <chrono>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/utils.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/gepetto/viewer.hpp"
#include "gepetto/viewer/corba/client.hh"

namespace idocp {

TrajectoryViewer::TrajectoryViewer(
    const std::string& description_pkg_serach_path, 
    const std::string& path_to_urdf)
  : description_pkg_serach_path_(description_pkg_serach_path),
    path_to_urdf_(path_to_urdf) {
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
  gui->setLightingMode("hpp-gui/floor", "OFF");

  for (const auto& q : q_traj) {
    viewer.display(q);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}

} // namespace idocp