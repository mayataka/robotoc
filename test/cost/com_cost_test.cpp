#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/com_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class CoMCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpulseCost(Robot& robot) const;

  double dt, t;
};


void CoMCostTest::testStageCost(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  cost->set_com_ref(com_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * com_weight.asDiagonal() * q_diff;
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, t, dt, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, t, dt, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += dt * J_3d.transpose() * com_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_3d.transpose() * com_weight.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


void CoMCostTest::testTerminalCost(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  cost->set_com_ref(com_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = 0.5 * q_diff.transpose() * comf_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, t, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * comf_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * comf_weight.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void CoMCostTest::testImpulseCost(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  cost->set_com_ref(com_ref);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = 0.5 * q_diff.transpose() * comi_weight.asDiagonal() * q_diff;
  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, t, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, t, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, t, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * comi_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * comi_weight.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


TEST_F(CoMCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(CoMCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}