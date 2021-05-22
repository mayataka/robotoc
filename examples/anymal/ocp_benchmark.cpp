#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/solver/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"

#include "idocp/utils/ocp_benchmarker.hpp"


int main () {
  // Create a robot with contacts.
  const int LF_foot = 12;
  const int LH_foot = 22;
  const int RF_foot = 32;
  const int RH_foot = 42;
  std::vector<int> contact_frames = {LF_foot, LH_foot, RF_foot, RH_foot}; // LF, LH, RF, RH
  const double baumgarte_time_step = 0.5 / 20;
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  idocp::Robot robot(path_to_urdf, idocp::BaseJointType::FloatingBase, 
                     contact_frames, baumgarte_time_step);

  // Create a cost function.
  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 0.4792, 0, 0, 0, 1, 
           -0.1,  0.7, -1.0, 
           -0.1, -0.7,  1.0, 
            0.1,  0.7, -1.0, 
            0.1, -0.7,  1.0;
  Eigen::VectorXd v_ref(robot.dimv());
  v_ref << 0, 0, 0, 0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0;
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_q_ref(q_ref);
  config_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    fr << 0, 0, 70;
    f_ref.push_back(fr);
  }
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(config_cost);
  cost->push_back(contact_cost);

  // Create inequality constraints.
  auto constraints = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone         = std::make_shared<idocp::FrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);

  // Create OCPSolver
  const double T = 0.5;
  const int N = 20;
  const int max_num_impulse_phase = 4;
  const int nthreads = 4;
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, 
                              max_num_impulse_phase, nthreads);

  // Initial time and initial state
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, 
       -0.1, -0.7,  1.0, 
        0.1,  0.7, -1.0, 
        0.1, -0.7,  1.0;
  Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  v << 0, 0, 0, 0, 0, 0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0;

  // Initialize OCPSolver
  auto contact_status = robot.createContactStatus();
  contact_status.activateContacts({0, 1, 2, 3});
  robot.updateFrameKinematics(q);
  robot.getContactPoints(contact_status);
  ocp_solver.setContactStatusUniformly(contact_status);
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  ocp_solver.initConstraints(t);

  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, 10, false);
  idocp::ocpbenchmarker::CPUTime(ocp_solver, t, q, v, 10000, false);

  // robot.printRobotModel();

  return 0;
}
