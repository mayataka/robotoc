#ifndef IDOCP_UTILS_SIMULATOR_HXX_
#define IDOCP_UTILS_SIMULATOR_HXX_

#include <chrono>
#include <assert.h>

namespace idocp {

namespace raisimwrapper {
inline void pino2rai(const Eigen::VectorXd& q_pinocchio, 
                     const Eigen::VectorXd& v_pinocchio, 
                     Eigen::VectorXd& q_raisim, Eigen::VectorXd& v_raisim) {
  q_raisim.coeffRef(0)  = q_pinocchio.coeff(0);
  q_raisim.coeffRef(1)  = q_pinocchio.coeff(1);
  q_raisim.coeffRef(2)  = q_pinocchio.coeff(2);
  q_raisim.coeffRef(3)  = q_pinocchio.coeff(6);
  q_raisim.coeffRef(4)  = q_pinocchio.coeff(3);
  q_raisim.coeffRef(5)  = q_pinocchio.coeff(4);
  q_raisim.coeffRef(6)  = q_pinocchio.coeff(5);
  q_raisim.coeffRef(7)  = q_pinocchio.coeff(7);
  q_raisim.coeffRef(8)  = q_pinocchio.coeff(8);
  q_raisim.coeffRef(9)  = q_pinocchio.coeff(9);
  q_raisim.coeffRef(10) = q_pinocchio.coeff(13);
  q_raisim.coeffRef(11) = q_pinocchio.coeff(14);
  q_raisim.coeffRef(12) = q_pinocchio.coeff(15);
  q_raisim.coeffRef(13) = q_pinocchio.coeff(10);
  q_raisim.coeffRef(14) = q_pinocchio.coeff(11);
  q_raisim.coeffRef(15) = q_pinocchio.coeff(12);
  q_raisim.coeffRef(16) = q_pinocchio.coeff(16);
  q_raisim.coeffRef(17) = q_pinocchio.coeff(17);
  q_raisim.coeffRef(18) = q_pinocchio.coeff(18);
  v_raisim.coeffRef(0)  = v_pinocchio.coeff(0);
  v_raisim.coeffRef(1)  = v_pinocchio.coeff(1);
  v_raisim.coeffRef(2)  = v_pinocchio.coeff(2);
  v_raisim.coeffRef(3)  = v_pinocchio.coeff(3);
  v_raisim.coeffRef(4)  = v_pinocchio.coeff(4);
  v_raisim.coeffRef(5)  = v_pinocchio.coeff(5);
  v_raisim.coeffRef(6)  = v_pinocchio.coeff(6);
  v_raisim.coeffRef(7)  = v_pinocchio.coeff(7);
  v_raisim.coeffRef(8)  = v_pinocchio.coeff(8);
  v_raisim.coeffRef(9)  = v_pinocchio.coeff(12);
  v_raisim.coeffRef(10) = v_pinocchio.coeff(13);
  v_raisim.coeffRef(11) = v_pinocchio.coeff(14);
  v_raisim.coeffRef(12) = v_pinocchio.coeff(9);
  v_raisim.coeffRef(13) = v_pinocchio.coeff(10);
  v_raisim.coeffRef(14) = v_pinocchio.coeff(11);
  v_raisim.coeffRef(15) = v_pinocchio.coeff(15);
  v_raisim.coeffRef(16) = v_pinocchio.coeff(16);
  v_raisim.coeffRef(17) = v_pinocchio.coeff(17);
}


inline void pino2rai(const Eigen::VectorXd& u_pino, Eigen::VectorXd& u_rai) {
  u_rai.template head<6>().setZero();
  u_rai.coeffRef(6)  = u_pino.coeff(6);
  u_rai.coeffRef(7)  = u_pino.coeff(7);
  u_rai.coeffRef(8)  = u_pino.coeff(8);
  u_rai.coeffRef(9)  = u_pino.coeff(12);
  u_rai.coeffRef(10) = u_pino.coeff(13);
  u_rai.coeffRef(11) = u_pino.coeff(14);
  u_rai.coeffRef(12) = u_pino.coeff(9);
  u_rai.coeffRef(13) = u_pino.coeff(10);
  u_rai.coeffRef(14) = u_pino.coeff(11);
  u_rai.coeffRef(15) = u_pino.coeff(15);
  u_rai.coeffRef(16) = u_pino.coeff(16);
  u_rai.coeffRef(17) = u_pino.coeff(17);
}


inline void rai2pino(const Eigen::VectorXd& q_raisim, 
                     const Eigen::VectorXd& v_raisim, 
                     Eigen::VectorXd& q_pinocchio,  
                     Eigen::VectorXd& v_pinocchio) {
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
  v_pinocchio.coeffRef(0)  = v_raisim.coeff(0);
  v_pinocchio.coeffRef(1)  = v_raisim.coeff(1);
  v_pinocchio.coeffRef(2)  = v_raisim.coeff(2);
  v_pinocchio.coeffRef(3)  = v_raisim.coeff(3);
  v_pinocchio.coeffRef(4)  = v_raisim.coeff(4);
  v_pinocchio.coeffRef(5)  = v_raisim.coeff(5);
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


inline void SetupCallback() {
  auto vis = raisim::OgreVis::get();
  /// light
  vis->getLight()->setDiffuseColour(1, 1, 1);
  vis->getLight()->setCastShadows(true);
  Ogre::Vector3 lightdir(-3,3,-0.5);
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
  // speed of camera motion in freelook mode
  vis->getCameraMan()->setTopSpeed(5);
  /// skybox
  Ogre::Quaternion quat;
  quat.FromAngleAxis(Ogre::Radian(0.5*M_PI_2), {1., 0, 0});
  vis->getSceneManager()->setSkyBox(true, "Examples/StormySkyBox", 500, true, quat);
}

} // namespace raisimwrapper


inline QuadrupedSimulator::QuadrupedSimulator(
    const std::string& path_to_urdf_for_raisim, 
    const std::string& save_dir_path, const std::string& save_file_name)
  : data_saver_(save_dir_path, save_file_name),
    raisim_world_(),
    raisim_robot_(raisim_world_.addArticulatedSystem(path_to_urdf_for_raisim, "")),
    raisim_ground_(raisim_world_.addGround()) {
}


template<typename OCPTypeDerived>
inline void QuadrupedSimulator::run(MPC<OCPTypeDerived>& mpc, 
                                    const double simulation_time_in_sec, 
                                    const double sampling_period_in_sec, 
                                    const double simulation_start_time_in_sec, 
                                    const Eigen::VectorXd& q_initial, 
                                    const Eigen::VectorXd& v_initial) {
  assert(simulation_time_in_sec > 0);
  assert(sampling_period_in_sec > 0);
  assert(simulation_start_time_in_sec >= 0);
  raisim_world_.setTimeStep(sampling_period_in_sec);
  raisim_world_.setERP(0.2, 0.2);

  Eigen::VectorXd q = q_initial;
  Eigen::VectorXd v = v_initial;
  Eigen::VectorXd u = Eigen::VectorXd::Zero(v_initial.size());
  Eigen::VectorXd q_raisim = Eigen::VectorXd::Zero(q_initial.size());
  Eigen::VectorXd v_raisim = Eigen::VectorXd::Zero(v_initial.size());
  Eigen::VectorXd u_raisim = Eigen::VectorXd::Zero(v_initial.size());

  mpc.getControlInput(u);
  raisimwrapper::pino2rai(q, v, q_raisim, v_raisim);
  raisimwrapper::pino2rai(u, u_raisim);

  std::chrono::system_clock::time_point start_clock, end_clock;
  double CPU_time_total_in_sec = 0;
  for (double t=0; t<simulation_time_in_sec; t+=sampling_period_in_sec) {
    raisimwrapper::pino2rai(q, v, q_raisim, v_raisim);
    raisimwrapper::pino2rai(u, u_raisim);
    raisim_robot_->setState(q_raisim, v_raisim);
    raisim_robot_->setGeneralizedForce(u_raisim);
    raisim_world_.integrate();
    data_saver_.save(q, v, u, mpc.computeKKTError(t, q, v));
    start_clock = std::chrono::system_clock::now();
    mpc.setContactPointByKinematics(q);
    mpc.updateSolution(t, q, v);
    mpc.getControlInput(u);
    end_clock = std::chrono::system_clock::now();
    const double CPU_time_in_sec
        = 1.0e-06 * std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    CPU_time_total_in_sec += CPU_time_in_sec;
    raisim_robot_->getState(q_raisim, v_raisim);
    raisimwrapper::rai2pino(q_raisim, v_raisim, q, v);
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
}


template<typename OCPTypeDerived>
inline void QuadrupedSimulator::viz(MPC<OCPTypeDerived>& mpc, 
                                    const double simulation_time_in_sec, 
                                    const double sampling_period_in_sec, 
                                    const double simulation_start_time_in_sec, 
                                    const Eigen::VectorXd& q_initial, 
                                    const Eigen::VectorXd& v_initial) {
  assert(simulation_time_in_sec > 0);
  assert(sampling_period_in_sec > 0);
  assert(simulation_start_time_in_sec >= 0);
  raisim_world_.setTimeStep(sampling_period_in_sec);
  raisim_world_.setERP(0.2, 0.2);

  auto vis = raisim::OgreVis::get();
  vis->setWorld(&raisim_world_);
  vis->setSetUpCallback(raisimwrapper::SetupCallback);
  vis->setAntiAliasing(2);
  vis->initApp();
  vis->setContactVisObjectSize(0.025, 0.01); 
  vis->setWindowSize(6000, 4000);
  vis->createGraphicalObject(raisim_ground_, 20, "floor", "checkerboard_green");
  vis->createGraphicalObject(raisim_robot_, "ANYmal");
  const double scaled_position = 0.6;
  vis->getCameraMan()->getCamera()->setPosition(-3.5*scaled_position, 
                                                -3*scaled_position, 
                                                1.5*scaled_position);
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(1.0, 0, 0), 
                                           Ogre::Radian(M_PI_2));
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(-1.0, -1.0, 0.), 
                                           Ogre::Radian(0.3));
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(0, -1.0, 0), 
                                           Ogre::Radian(0.5));
  vis->getCameraMan()->getCamera()->rotate(Ogre::Vector3(0, 0, -1.0), 
                                           Ogre::Radian(0.2));

  Eigen::VectorXd q = q_initial;
  Eigen::VectorXd v = v_initial;
  Eigen::VectorXd u = Eigen::VectorXd::Zero(v_initial.size());
  Eigen::VectorXd q_raisim = Eigen::VectorXd::Zero(q_initial.size());
  Eigen::VectorXd v_raisim = Eigen::VectorXd::Zero(v_initial.size());
  Eigen::VectorXd u_raisim = Eigen::VectorXd::Zero(v_initial.size());

  mpc.getControlInput(u);
  raisimwrapper::pino2rai(q, v, q_raisim, v_raisim);
  raisimwrapper::pino2rai(u, u_raisim);

  std::chrono::system_clock::time_point start_clock, end_clock;
  double CPU_time_total_in_sec = 0;
  for (double t=0; t<simulation_time_in_sec; t+=sampling_period_in_sec) {
    raisimwrapper::pino2rai(q, v, q_raisim, v_raisim);
    raisimwrapper::pino2rai(u, u_raisim);
    raisim_robot_->setState(q_raisim, v_raisim);
    raisim_robot_->setGeneralizedForce(u_raisim);
    raisim_world_.integrate();
    vis->renderOneFrame();
    data_saver_.save(q, v, u, mpc.computeKKTError(t, q, v));
    start_clock = std::chrono::system_clock::now();
    mpc.setContactPointByKinematics(q);
    mpc.updateSolution(t, q, v);
    mpc.getControlInput(u);
    end_clock = std::chrono::system_clock::now();
    const double CPU_time_in_sec
        = 1.0e-06 * std::chrono::duration_cast<std::chrono::microseconds>(
            end_clock-start_clock).count();
    CPU_time_total_in_sec += CPU_time_in_sec;
    raisim_robot_->getState(q_raisim, v_raisim);
    raisimwrapper::rai2pino(q_raisim, v_raisim, q, v);
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

  vis->closeApp();
}

} // namespace idocp

#endif // IDOCP_UTILS_SIMULATOR_HXX_ 