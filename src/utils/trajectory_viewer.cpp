#include "idocp/utils/trajectory_viewer.hpp"
#include "idocp/utils/helper.hpp"
#include "idocp/utils/csv_to_eigen.hpp"

#include "raisim/World.hpp"
#include "raisim/OgreVis.hpp"


namespace idocp {

TrajectoryViewer::TrajectoryViewer(
    const std::string& path_to_raisim_activation_key,
    const std::string& path_to_urdf_for_raisim) 
  : path_to_raisim_activation_key_(path_to_raisim_activation_key),
    path_to_urdf_for_raisim_(path_to_urdf_for_raisim) {
}


void TrajectoryViewer::display(const std::string& path_to_data, 
                               const double sampling_period_in_sec, 
                               const bool recording, 
                               const std::string& path_to_save_file) {
  raisim::World::setActivationKey(path_to_raisim_activation_key_);
  raisim::World raisim_world;
  auto raisim_robot = raisim_world.addArticulatedSystem(
      helper::GetCurrentWorkingDirectory()+path_to_urdf_for_raisim_, "");
  auto raisim_ground = raisim_world.addGround();
  raisim_world.setTimeStep(sampling_period_in_sec);
  auto vis = raisim::OgreVis::get();
  vis->setWorld(&raisim_world);
  vis->setSetUpCallback(TrajectoryViewer::RaiSimOgreCallback);
  vis->setAntiAliasing(2);
  vis->setWindowSize(500, 400);
  vis->setContactVisObjectSize(0.025, 0.01); 
  vis->initApp();
  vis->getCameraMan()->getCamera()->setPosition(2.5, 2.5, 1);
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(0.0, 0.0, 1.0), 
                                           Ogre::Radian(3*M_PI_4));
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(1.0, 0.0, 0.0), 
                                           Ogre::Radian(M_PI_2*9/10));
  vis->createGraphicalObject(raisim_ground, 20, "floor", "checkerboard_green");
  vis->createGraphicalObject(raisim_robot, "ANYmal");
  const int dimq = 19;
  CSVToEigen csv_to_eigen(path_to_data, dimq);
  csv_to_eigen.readCSVLine();
  raisim_robot->setGeneralizedCoordinate(helper::pino2rai_q(csv_to_eigen.get()));
  vis->renderOneFrame();
  if (recording) {
    vis->startRecordingVideo(path_to_save_file+".mp4");
  }
  vis->renderOneFrame();
  while (csv_to_eigen.readCSVLine()) {
    raisim_robot->setGeneralizedCoordinate(helper::pino2rai_q(csv_to_eigen.get()));
    vis->renderOneFrame();
  }
  if (recording) {
    vis->stopRecordingVideoAndSave();
  }
  vis->closeApp();
}


void TrajectoryViewer::RaiSimOgreCallback() {
  auto vis = raisim::OgreVis::get();
  /// light
  vis->getLight()->setDiffuseColour(1, 1, 1);
  vis->getLight()->setCastShadows(true);
  Ogre::Vector3 lightdir(-3, -3, -0.5);
  lightdir.normalise();
  vis->getLightNode()->setDirection({lightdir});
  vis->setCameraSpeed(300);
  /// load  textures
  vis->addResourceDirectory(vis->getResourceDir() + "/material/checkerboard");
  vis->loadMaterialFile("checkerboard.material");
  /// shdow setting
  vis->getSceneManager()->setShadowTechnique(Ogre::SHADOWTYPE_TEXTURE_ADDITIVE);
  vis->getSceneManager()->setShadowTextureSettings(2048, 3);
  /// scale related settings!! Please adapt it depending on your map size
  /// beyond this distance, shadow disappears
  vis->getSceneManager()->setShadowFarDistance(10);
  // size of contact points and contact forces
  // vis->setContactVisObjectSize(0.03, 0.6);
  vis->setContactVisObjectSize(0.03, 6);
  // speed of camera motion in free look mode
  vis->getCameraMan()->setTopSpeed(5);
  /// skybox
  Ogre::Quaternion quat;
  quat.FromAngleAxis(Ogre::Radian(M_PI_2), {1., 0, 0});
  vis->getSceneManager()->setSkyBox(true, "Examples/StormySkyBox", 500, true, quat);
  // vis->getSceneManager()->setSkyBox(true, "white", 500, true, quat);
}

} // namespace idocp