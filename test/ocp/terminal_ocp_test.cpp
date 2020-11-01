#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
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

  static void testLinearizeOCPAndBackwardRiccatiRecursion(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testTerminalCost(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testComputeCondensedPrimalDirection(
      const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testComputeCondensedDualDirection(
      const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testUpdatePrimal(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testUpdateDual(
      const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  static void testComputeKKTResidual(
      Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints);

  std::string fixed_base_urdf, floating_base_urdf;
};


std::shared_ptr<CostFunction> TerminalOCPTest::createCost(const Robot& robot) {
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimu()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_vf_weight(vf_weight);
  const int task_frame = 10;
  auto task_space_3d_cost = std::make_shared<TaskSpace3DCost >(robot, task_frame);
  const Eigen::Vector3d q_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d q_3d_ref = Eigen::Vector3d::Random();
  task_space_3d_cost->set_q_3d_weight(q_3d_weight);
  task_space_3d_cost->set_qf_3d_weight(qf_3d_weight);
  task_space_3d_cost->set_q_3d_ref(q_3d_ref);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(joint_cost);
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


void TerminalOCPTest::testLinearizeOCPAndBackwardRiccatiRecursion(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  ocp.linearizeOCP(robot, t, s);
  RiccatiSolution riccati(robot);
  robot.updateKinematics(s.q, s.v);
  ocp.backwardRiccatiRecursion(riccati);
  KKTResidual kkt_residual(robot);  
  KKTMatrix kkt_matrix(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  cost->computeTerminalCostHessian(robot, cost_data, t, s, kkt_matrix);
  EXPECT_TRUE(riccati.sq.isApprox(-1*kkt_residual.lq()+s.lmd));
  EXPECT_TRUE(riccati.sv.isApprox(-1*kkt_residual.lv()+s.gmm));
  EXPECT_TRUE(riccati.Pqq.isApprox(kkt_matrix.Qqq()));
  EXPECT_TRUE(riccati.Pqv.isApprox(kkt_matrix.Qqv()));
  EXPECT_TRUE(riccati.Pvq.isApprox(kkt_matrix.Qvq()));
  EXPECT_TRUE(riccati.Pvv.isApprox(kkt_matrix.Qvv()));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati.Pvq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  double KKT_ref = 0;
  KKT_ref += (-1*kkt_residual.lq()+s.lmd).squaredNorm();
  KKT_ref += (-1*kkt_residual.lv()+s.gmm).squaredNorm();
  EXPECT_DOUBLE_EQ(KKT_ref, ocp.squaredNormKKTResidual());
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
  const double terminal_cost_ref = cost->phi(robot, cost_data, t, s);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
}


void TerminalOCPTest::testComputeCondensedPrimalDirection(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  RiccatiSolution riccati(robot);
  riccati.Pqq.setRandom();
  riccati.Pqv.setRandom();
  riccati.Pvq.setRandom();
  riccati.Pvv.setRandom();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  SplitDirection d(robot);
  d.dq().setRandom();
  d.dv().setRandom();
  SplitDirection d_ref = d;
  ocp.computeCondensedPrimalDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d_ref.dq() + riccati.Pqv * d_ref.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d_ref.dq() + riccati.Pvv * d_ref.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void TerminalOCPTest::testComputeCondensedDualDirection(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  TerminalOCP ocp(robot, cost, constraints);
  const SplitDirection d = SplitDirection::Random(robot);
  ocp.computeCondensedDualDirection(d);
}


void TerminalOCPTest::testUpdatePrimal(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_ref = s;
  const SplitDirection d = SplitDirection::Random(robot);
  const double step_size = 0.3;
  TerminalOCP ocp(robot, cost, constraints);
  ocp.updatePrimal(robot, step_size, d, s);
  s_ref.lmd += step_size * d.dlmd();
  s_ref.gmm += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s_ref.q);
  s_ref.v += step_size * d.dv();
  EXPECT_TRUE(s.isApprox(s_ref));
}


void TerminalOCPTest::testUpdateDual(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const SplitDirection d = SplitDirection::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  const double step_size = 0.3;
  ocp.updateDual(step_size);
}


void TerminalOCPTest::testComputeKKTResidual(
    Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  TerminalOCP ocp(robot, cost, constraints);
  ocp.computeKKTResidual(robot, t, s);
  const double KKT = ocp.squaredNormKKTResidual();
  robot.updateKinematics(s.q, s.v);
  KKTResidual kkt_residual(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  double KKT_ref = 0;
  KKT_ref += (-1*kkt_residual.lq()+s.lmd).squaredNorm();
  KKT_ref += (-1*kkt_residual.lv()+s.gmm).squaredNorm();
  EXPECT_DOUBLE_EQ(KKT, KKT_ref);
}


TEST_F(TerminalOCPTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCPAndBackwardRiccatiRecursion(robot, cost, constraints);
  testTerminalCost(robot, cost, constraints);
  testComputeCondensedPrimalDirection(robot, cost, constraints);
  testComputeCondensedDualDirection(robot, cost, constraints);
  testUpdatePrimal(robot, cost, constraints);
  testUpdateDual(robot, cost, constraints);
  testComputeKKTResidual(robot, cost, constraints);
}


TEST_F(TerminalOCPTest, floatingBase) {
  Robot robot(floating_base_urdf);
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCPAndBackwardRiccatiRecursion(robot, cost, constraints);
  testTerminalCost(robot, cost, constraints);
  testComputeCondensedPrimalDirection(robot, cost, constraints);
  testComputeCondensedDualDirection(robot, cost, constraints);
  testUpdatePrimal(robot, cost, constraints);
  testUpdateDual(robot, cost, constraints);
  testComputeKKTResidual(robot, cost, constraints);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}