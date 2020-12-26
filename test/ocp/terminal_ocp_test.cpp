#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"


namespace idocp {

class TerminalOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);

  static std::shared_ptr<Constraints> createConstraints(const Robot& robot);

  static void testLinearizeOCP(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testTerminalCost(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testComputeKKTResidual(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  std::string fixed_base_urdf, floating_base_urdf;
};


std::shared_ptr<CostFunction> TerminalOCPTest::createCost(const Robot& robot) {
  auto config_cost = std::make_shared<ConfigurationSpaceCost >(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimu()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_ref(v_ref);
  config_cost->set_a_weight(a_weight);
  config_cost->set_u_weight(u_weight);
  config_cost->set_u_ref(u_ref);
  config_cost->set_qf_weight(qf_weight);
  config_cost->set_vf_weight(vf_weight);
  const int task_frame = 10;
  auto task_space_3d_cost = std::make_shared<TaskSpace3DCost >(robot, task_frame);
  const Eigen::Vector3d q_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d q_3d_ref = Eigen::Vector3d::Random();
  task_space_3d_cost->set_q_3d_weight(q_3d_weight);
  task_space_3d_cost->set_qf_3d_weight(qf_3d_weight);
  task_space_3d_cost->set_q_3d_ref(q_3d_ref);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(config_cost);
  cost->push_back(task_space_3d_cost);
  return cost;
}


std::shared_ptr<Constraints> TerminalOCPTest::createConstraints(const Robot& robot) {
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto constraints = std::make_shared<Constraints>();
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  return constraints;
}


void TerminalOCPTest::testLinearizeOCP(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTMatrix kkt_matrix(robot);  
  SplitKKTResidual kkt_residual(robot);  
  ocp.linearizeOCP(robot, t, s, kkt_matrix, kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTMatrix kkt_matrix_ref(robot);  
  SplitKKTResidual kkt_residual_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  kkt_residual_ref.lq() -= s.lmd;
  kkt_residual_ref.lv() -= s.gmm;
  cost->computeTerminalCostHessian(robot, cost_data, t, s, kkt_matrix_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void TerminalOCPTest::testTerminalCost(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  const double terminal_cost = ocp.terminalCost(robot, t, s);
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost_ref = cost->computeTerminalCost(robot, cost_data, t, s);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
}


void TerminalOCPTest::testComputeKKTResidual(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  ocp.computeKKTResidual(robot, t, s, kkt_residual);
  const double KKT = ocp.squaredNormKKTResidual(kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTResidual kkt_residual_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  kkt_residual_ref.lq() -= s.lmd;
  kkt_residual_ref.lv() -= s.gmm;
  double KKT_ref = 0;
  KKT_ref += kkt_residual_ref.lq().squaredNorm();
  KKT_ref += kkt_residual_ref.lv().squaredNorm();
  EXPECT_DOUBLE_EQ(KKT, KKT_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(TerminalOCPTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, cost, constraints);
  testTerminalCost(robot, cost, constraints);
  testComputeKKTResidual(robot, cost, constraints);
}


TEST_F(TerminalOCPTest, floatingBase) {
  Robot robot(floating_base_urdf);
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, cost, constraints);
  testTerminalCost(robot, cost, constraints);
  testComputeKKTResidual(robot, cost, constraints);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}