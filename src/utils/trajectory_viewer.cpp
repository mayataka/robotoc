#include "idocp/utils/trajectory_viewer.hpp"

#include <thread>
#include <chrono>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/utils.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/gepetto/viewer.hpp"
#include "gepetto/viewer/corba/client.hh"

namespace idocp {

TrajectoryViewer::TrajectoryViewer(const std::string& path_to_description_pkg, 
                                   const std::string& path_to_urdf)
  : path_to_description_pkg_(path_to_description_pkg),
    path_to_urdf_(path_to_urdf) {
}


void TrajectoryViewer::display(const std::vector<Eigen::VectorXd>& q_traj, 
                               const double sampling_period_in_sec, 
                               const bool recording, 
                               const std::string& path_to_save_file) {
  pinocchio::Model model;
  pinocchio::urdf::buildModel(path_to_urdf_, model);
  pinocchio::GeometryModel vmodel;
  pinocchio::urdf::buildGeom(model, path_to_urdf_, pinocchio::VISUAL, vmodel, 
                             path_to_description_pkg_);
  
  pinocchio::gepetto::Viewer viewer(model, &vmodel, NULL);
  viewer.initViewer("idocp::TrajectoryViewer", true);
  viewer.loadViewerModel("world");

  auto gui = gepetto::viewer::corba::gui();
  gui->setColor("hpp-gui/floor", gepetto::viewer::corba::grey);
  gui->setLightingMode("hpp-gui/floor", "OFF");

  for (const auto& q : q_traj) {
    viewer.display(q);
    std::this_thread::sleep_for(std::chrono::duration<double>(sampling_period_in_sec));
  }
}

} // namespace idocp