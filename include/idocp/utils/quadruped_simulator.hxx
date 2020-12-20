#ifndef IDOCP_UTILS_QUADRUPED_SIMULATOR_HXX_
#define IDOCP_UTILS_QUADRUPED_SIMULATOR_HXX_ 

#include "idocp/utils/quadruped_simulator.hpp"

#include <chrono>
#include <stdexcept>

namespace idocp {

template <typename OCPSolverType>
inline QuadrupedSimulator<OCPSolverType>::QuadrupedSimulator(
    const std::string& path_to_raisim_activation_key,
    const std::string& path_to_urdf_for_raisim, 
    const std::string& save_dir_path, const std::string& save_file_name)
  : path_to_raisim_activation_key_(path_to_raisim_activation_key),
    path_to_urdf_for_raisim_(path_to_urdf_for_raisim),
    data_saver_(save_dir_path, save_file_name),
    save_dir_path_(save_dir_path), 
    save_file_name_(save_file_name) {
}


template <typename OCPSolverType>
template <typename MPCCallbackType>
inline void QuadrupedSimulator<OCPSolverType>::run(
    MPC<OCPSolverType>& mpc, const double simulation_time_in_sec, 
    const double sampling_period_in_sec, 
    const double simulation_start_time_in_sec, 
    const Eigen::VectorXd& q_initial, const Eigen::VectorXd& v_initial, 
    const bool visualization, const bool recording) {
  try {
    if (simulation_time_in_sec <= 0) {
      throw std::out_of_range(
          "Invalid argument: simulation_time_in_sec must be positive!");
    }
    if (sampling_period_in_sec <= 0) {
      throw std::out_of_range(
          "Invalid argument: sampling_period_in_sec must be positive!");
    }
    if (simulation_start_time_in_sec < 0) {
      throw std::out_of_range(
          "Invalid argument: simulation_start_time_in_sec must be non negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  raisim::World::setActivationKey(path_to_raisim_activation_key_);
  raisim::World raisim_world;
  auto raisim_robot = raisim_world.addArticulatedSystem(path_to_urdf_for_raisim_, "");
  auto raisim_ground = raisim_world.addGround();
  raisim_world.setTimeStep(sampling_period_in_sec);
  raisim_world.setERP(0.2, 0.2);
  // raisim_world.setDefaultMaterial(1000, 0, 0);
  auto vis = raisim::OgreVis::get();
  if (visualization) {
    raisimadapter::SetupRaisimOgre(raisim_world);
    vis->createGraphicalObject(raisim_ground, 20, "floor", "checkerboard_green");
    vis->createGraphicalObject(raisim_robot, "ANYmal");
  }
  Robot robot = mpc.getSolverHandle()->createRobot();
  Eigen::VectorXd q = q_initial;
  Eigen::VectorXd v = v_initial;
  Eigen::VectorXd u = Eigen::VectorXd::Zero(robot.dimu());
  Eigen::VectorXd q_raisim = Eigen::VectorXd::Zero(robot.dimq());
  Eigen::VectorXd v_raisim = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd u_raisim = Eigen::VectorXd::Zero(robot.dimv());
  mpc.getControlInput(u);
  raisimadapter::pino2rai(robot, q, v, q_raisim, v_raisim);
  raisim_robot->setState(q_raisim, v_raisim);
  raisim_robot->setGeneralizedForce(u_raisim);
  if (visualization && recording) {
    vis->startRecordingVideo(save_dir_path_+'/'+save_file_name_+".mp4");
  }
  std::chrono::system_clock::time_point start_clock, end_clock;
  double CPU_time_total_in_sec = 0;
  MPCCallbackType mpc_callback(robot);
  mpc_callback.init(simulation_start_time_in_sec, q, v, mpc);
  for (double t=simulation_start_time_in_sec; t<simulation_time_in_sec; t+=sampling_period_in_sec) {
    raisim_robot->setGeneralizedForce(u_raisim);
    raisim_world.integrate();
    if (visualization) {
      vis->renderOneFrame();
    }
    raisim_robot->getState(q_raisim, v_raisim);
    raisimadapter::rai2pino(robot, q_raisim, v_raisim, q, v);
    mpc.computeKKTResidual(t, q, v);
    data_saver_.save(q, v, u, mpc.KKTError());
    start_clock = std::chrono::system_clock::now();
    mpc_callback.callback(t, q, v, mpc);
    mpc.updateSolution(t, q, v);
    mpc.getControlInput(u);
    raisimadapter::pino2rai(u, u_raisim);
    end_clock = std::chrono::system_clock::now();
    const double CPU_time_in_sec
        = 1.0e-06 * std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    CPU_time_total_in_sec += CPU_time_in_sec;
  }
  const int num_simulation_update 
      = (int)(simulation_time_in_sec/sampling_period_in_sec);
  const double CPU_time_per_update_in_sec 
      = CPU_time_total_in_sec / num_simulation_update;
  data_saver_.saveConditions(simulation_time_in_sec, 
                              1000*sampling_period_in_sec,
                              1000*CPU_time_per_update_in_sec);
  std::cout << "simulation time: " << simulation_time_in_sec << "[s]" 
            << std::endl;
  std::cout << "sampling time: " << 1000 * sampling_period_in_sec << "[ms]" 
            << std::endl;
  std::cout << "CPU time per update: " << 1000 * CPU_time_per_update_in_sec 
            << "[ms]" << std::endl;
  if (visualization && recording) {
    vis->stopRecordingVideoAndSave();
  }
  vis->closeApp();
}


inline void raisimadapter::SetupRaisimOgre(raisim::World& world) {
  auto vis = raisim::OgreVis::get();
  vis->setWorld(&world);
  vis->setSetUpCallback(raisimadapter::RaiSimOgreCallback);
  vis->setAntiAliasing(2);
  vis->setWindowSize(500, 400);
  vis->setContactVisObjectSize(0.025, 0.01); 
  vis->initApp();
  vis->getCameraMan()->getCamera()->setPosition(2.5, 2.5, 1);
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(0.0, 0.0, 1.0), 
                                           Ogre::Radian(3*M_PI_4));
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(1.0, 0.0, 0.0), 
                                           Ogre::Radian(M_PI_2*9/10));
}


inline void raisimadapter::RaiSimOgreCallback() {
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
  vis->setContactVisObjectSize(0.03, 0.6);
  // speed of camera motion in free look mode
  vis->getCameraMan()->setTopSpeed(5);
  /// skybox
  Ogre::Quaternion quat;
  quat.FromAngleAxis(Ogre::Radian(M_PI_2), {1., 0, 0});
  vis->getSceneManager()->setSkyBox(true, "Examples/StormySkyBox", 500, true, quat);
  // vis->getSceneManager()->setSkyBox(true, "white", 500, true, quat);
}


inline void raisimadapter::pino2rai(Robot& robot, 
                                    const Eigen::VectorXd& q_pinocchio, 
                                    const Eigen::VectorXd& v_pinocchio, 
                                    Eigen::VectorXd& q_raisim, 
                                    Eigen::VectorXd& v_raisim) {
  assert(q_pinocchio.size() == 19);
  assert(v_pinocchio.size() == 18);
  assert(q_raisim.size() == 19);
  assert(v_raisim.size() == 18);
  q_raisim.coeffRef(0)  = q_pinocchio.coeff(0);
  q_raisim.coeffRef(1)  = q_pinocchio.coeff(1);
  q_raisim.coeffRef(2)  = q_pinocchio.coeff(2);
  // Change quaternion order
  q_raisim.coeffRef(3)  = q_pinocchio.coeff(6);
  q_raisim.coeffRef(4)  = q_pinocchio.coeff(3);
  q_raisim.coeffRef(5)  = q_pinocchio.coeff(4);
  q_raisim.coeffRef(6)  = q_pinocchio.coeff(5);
  q_raisim.coeffRef(7)  = q_pinocchio.coeff(7);
  q_raisim.coeffRef(8)  = q_pinocchio.coeff(8);
  q_raisim.coeffRef(9)  = q_pinocchio.coeff(9);
  // Change leg order
  q_raisim.coeffRef(10) = q_pinocchio.coeff(13);
  q_raisim.coeffRef(11) = q_pinocchio.coeff(14);
  q_raisim.coeffRef(12) = q_pinocchio.coeff(15);
  // Change leg order
  q_raisim.coeffRef(13) = q_pinocchio.coeff(10);
  q_raisim.coeffRef(14) = q_pinocchio.coeff(11);
  q_raisim.coeffRef(15) = q_pinocchio.coeff(12);
  q_raisim.coeffRef(16) = q_pinocchio.coeff(16);
  q_raisim.coeffRef(17) = q_pinocchio.coeff(17);
  q_raisim.coeffRef(18) = q_pinocchio.coeff(18);
  // Change velocity's reference frame from LOCAL to WORLD.
  robot.updateFrameKinematics(q_pinocchio);
  const Eigen::Matrix3d R = robot.frameRotation(4);
  v_raisim.head<3>()     = R * v_pinocchio.head<3>();
  v_raisim.segment<3>(3) = R * v_pinocchio.segment<3>(3);
  v_raisim.coeffRef(6)  = v_pinocchio.coeff(6);
  v_raisim.coeffRef(7)  = v_pinocchio.coeff(7);
  v_raisim.coeffRef(8)  = v_pinocchio.coeff(8);
  // Change leg order
  v_raisim.coeffRef(9)  = v_pinocchio.coeff(12);
  v_raisim.coeffRef(10) = v_pinocchio.coeff(13);
  v_raisim.coeffRef(11) = v_pinocchio.coeff(14);
  // Change leg order
  v_raisim.coeffRef(12) = v_pinocchio.coeff(9);
  v_raisim.coeffRef(13) = v_pinocchio.coeff(10);
  v_raisim.coeffRef(14) = v_pinocchio.coeff(11);
  v_raisim.coeffRef(15) = v_pinocchio.coeff(15);
  v_raisim.coeffRef(16) = v_pinocchio.coeff(16);
  v_raisim.coeffRef(17) = v_pinocchio.coeff(17);
}


inline void raisimadapter::pino2rai(const Eigen::VectorXd& u_pinocchio, 
                                    Eigen::VectorXd& u_raisim) {
  assert(u_pinocchio.size() == 12);
  assert(u_raisim.size() == 18);
  u_raisim.template head<6>().setZero();
  u_raisim.coeffRef(6)  = u_pinocchio.coeff(0);
  u_raisim.coeffRef(7)  = u_pinocchio.coeff(1);
  u_raisim.coeffRef(8)  = u_pinocchio.coeff(2);
  // Change leg order
  u_raisim.coeffRef(9)  = u_pinocchio.coeff(6);
  u_raisim.coeffRef(10) = u_pinocchio.coeff(7);
  u_raisim.coeffRef(11) = u_pinocchio.coeff(8);
  // Change leg order
  u_raisim.coeffRef(12) = u_pinocchio.coeff(3);
  u_raisim.coeffRef(13) = u_pinocchio.coeff(4);
  u_raisim.coeffRef(14) = u_pinocchio.coeff(5);
  // Change leg order
  u_raisim.coeffRef(15) = u_pinocchio.coeff(9);
  u_raisim.coeffRef(16) = u_pinocchio.coeff(10);
  u_raisim.coeffRef(17) = u_pinocchio.coeff(11);
}


inline void raisimadapter::rai2pino(Robot& robot, 
                                    const Eigen::VectorXd& q_raisim, 
                                    const Eigen::VectorXd& v_raisim, 
                                    Eigen::VectorXd& q_pinocchio, 
                                    Eigen::VectorXd& v_pinocchio) {
  assert(q_raisim.size() == 19);
  assert(v_raisim.size() == 18);
  assert(q_pinocchio.size() == 19);
  assert(v_pinocchio.size() == 18);
  q_pinocchio.coeffRef(0)  = q_raisim.coeff(0);
  q_pinocchio.coeffRef(1)  = q_raisim.coeff(1);
  q_pinocchio.coeffRef(2)  = q_raisim.coeff(2);

  q_pinocchio.coeffRef(6)  = q_raisim.coeff(3);
  q_pinocchio.coeffRef(3)  = q_raisim.coeff(4);
  q_pinocchio.coeffRef(4)  = q_raisim.coeff(5);
  q_pinocchio.coeffRef(5)  = q_raisim.coeff(6);

  q_pinocchio.coeffRef(7)  = q_raisim.coeff(7);
  q_pinocchio.coeffRef(8)  = q_raisim.coeff(8);
  q_pinocchio.coeffRef(9)  = q_raisim.coeff(9);

  q_pinocchio.coeffRef(13) = q_raisim.coeff(10);
  q_pinocchio.coeffRef(14) = q_raisim.coeff(11);
  q_pinocchio.coeffRef(15) = q_raisim.coeff(12);

  q_pinocchio.coeffRef(10) = q_raisim.coeff(13);
  q_pinocchio.coeffRef(11) = q_raisim.coeff(14);
  q_pinocchio.coeffRef(12) = q_raisim.coeff(15);

  q_pinocchio.coeffRef(16) = q_raisim.coeff(16);
  q_pinocchio.coeffRef(17) = q_raisim.coeff(17);
  q_pinocchio.coeffRef(18) = q_raisim.coeff(18);

  // Change velocity's reference frame from WORLD to LOCAL.
  robot.updateFrameKinematics(q_pinocchio);
  const Eigen::Matrix3d R = robot.frameRotation(4);
  v_pinocchio.head<3>()     = R.transpose() * v_raisim.head<3>();
  v_pinocchio.segment<3>(3) = R.transpose() * v_raisim.segment<3>(3);

  v_pinocchio.coeffRef(6)  = v_raisim.coeff(6);
  v_pinocchio.coeffRef(7)  = v_raisim.coeff(7);
  v_pinocchio.coeffRef(8)  = v_raisim.coeff(8);

  v_pinocchio.coeffRef(9)  = v_raisim.coeff(12);
  v_pinocchio.coeffRef(10) = v_raisim.coeff(13);
  v_pinocchio.coeffRef(11) = v_raisim.coeff(14);

  v_pinocchio.coeffRef(12) = v_raisim.coeff(9);
  v_pinocchio.coeffRef(13) = v_raisim.coeff(10);
  v_pinocchio.coeffRef(14) = v_raisim.coeff(11);

  v_pinocchio.coeffRef(15) = v_raisim.coeff(15);
  v_pinocchio.coeffRef(16) = v_raisim.coeff(16);
  v_pinocchio.coeffRef(17) = v_raisim.coeff(17);
}

} // namespace idocp

#endif // IDOCP_UTILS_QUADRUPED_SIMULATOR_HXX_ 