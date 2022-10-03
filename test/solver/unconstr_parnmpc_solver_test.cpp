#include <vector>

#include <gtest/gtest.h>

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/solver/unconstr_parnmpc_solver.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class UnconstrParNMPCSolverTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(UnconstrParNMPCSolverTest, test) {
  auto robot = testhelper::CreateRobotManipulator();

  // // Change the limits from the default parameters.
  // robot.setJointEffortLimit(Eigen::VectorXd::Constant(robot.dimu(), 50));
  // robot.setJointVelocityLimit(Eigen::VectorXd::Constant(robot.dimv(), M_PI_2));

  // Create a cost function.
  auto cost = std::make_shared<robotoc::CostFunction>();
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_ref(Eigen::VectorXd::Zero(robot.dimq()));
  q_ref << 0, M_PI_2, 0, M_PI_2, 0, M_PI_2, 0;
  config_cost->set_q_ref(q_ref);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_q_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(config_cost);

  // Create joint constraints.
  const double barrier_param = 1.0e-03;
  const double fraction_to_boundary_rule = 0.995;
  auto constraints = std::make_shared<robotoc::Constraints>(barrier_param, fraction_to_boundary_rule);
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);

  // Create the ParNMPC solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  robotoc::OCP ocp(robot, cost, constraints, T, N);
  auto solver_options = robotoc::SolverOptions();
  const int nthreads = 4;
  robotoc::UnconstrParNMPCSolver ocp_solver(ocp, solver_options, nthreads);

  // Initial time and initial state
  const double t = 0;
  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << M_PI_2, 0, M_PI_2, 0, M_PI_2, 0, M_PI_2;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());

  // Solves the ParNMPC.
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  ocp_solver.initConstraints();
  ocp_solver.solve(t, q, v);
  const auto result = ocp_solver.getSolverStatistics();
  EXPECT_TRUE(result.convergence);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}