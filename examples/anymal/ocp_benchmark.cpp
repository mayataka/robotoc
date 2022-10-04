#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/local_contact_force_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/solver/solver_options.hpp"

#include "robotoc/utils/ocp_benchmarker.hpp"


int main () {
  // Create a robot with contacts.
  robotoc::RobotModelInfo model_info;
  model_info.urdf_path = "../anymal_b_simple_description/urdf/anymal.urdf";
  model_info.base_joint_type = robotoc::BaseJointType::FloatingBase;
  const double baumgarte_time_step = 0.025;
  model_info.point_contacts = {robotoc::ContactModelInfo("LF_FOOT", baumgarte_time_step),
                               robotoc::ContactModelInfo("LH_FOOT", baumgarte_time_step),
                               robotoc::ContactModelInfo("RF_FOOT", baumgarte_time_step),
                               robotoc::ContactModelInfo("RH_FOOT", baumgarte_time_step)};
  robotoc::Robot robot(model_info);

  // Create a cost function.
  auto cost = std::make_shared<robotoc::CostFunction>();
  Eigen::VectorXd q_standing(robot.dimq());
  q_standing << 0, 0, 0.4792, 0, 0, 0, 1, 
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
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_q_ref(q_standing);
  config_cost->set_q_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(config_cost);
  auto local_contact_force_cost = std::make_shared<robotoc::LocalContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<model_info.point_contacts.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    fr << 0, 0, 70;
    f_ref.push_back(fr);
  }
  local_contact_force_cost->set_f_weight(f_weight);
  local_contact_force_cost->set_f_ref(f_ref);
  cost->push_back(local_contact_force_cost);

  // Create inequality constraints.
  auto constraints = std::make_shared<robotoc::Constraints>();
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  auto friction_cone        = std::make_shared<robotoc::FrictionCone>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);

  // Create the contact sequence
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot);
  const double mu = 0.7;

  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts({0, 1, 2, 3});
  robot.updateFrameKinematics(q_standing);
  const std::vector<Eigen::Vector3d> contact_positions = {robot.framePosition("LF_FOOT"), 
                                                          robot.framePosition("LH_FOOT"),
                                                          robot.framePosition("RF_FOOT"),
                                                          robot.framePosition("RH_FOOT")};
  const std::vector<double> friction_coefficients = {mu, mu, mu, mu};
  contact_status_standing.setContactPlacements(contact_positions);
  contact_status_standing.setFrictionCoefficients(friction_coefficients);
  contact_sequence->init(contact_status_standing);

  // Create OCPSolver
  const double T = 0.5;
  const int N = 20;
  robotoc::OCP ocp(robot, cost, constraints, contact_sequence, T, N);
  auto solver_options = robotoc::SolverOptions();
  solver_options.nthreads = 4;
  robotoc::OCPSolver ocp_solver(ocp, solver_options);

  // Initial time and initial state
  const double t = 0;
  const Eigen::VectorXd q = q_standing;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());

  ocp_solver.discretize(t);
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  ocp_solver.initConstraints();
  ocp_solver.solve(t, q, v);
  std::cout << ocp_solver.getSolverStatistics() << std::endl;

  const int num_iteration = 10000;
  robotoc::benchmark::CPUTime(ocp_solver, t, q, v, num_iteration);

  // std::cout << robot << std::endl;

  return 0;
}
