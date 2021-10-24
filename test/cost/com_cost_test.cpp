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
  const Eigen::Vector3d q_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d CoM_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_weight(q_weight);
  cost->set_qf_weight(qf_weight);
  cost->set_qi_weight(qi_weight);
  cost->set_CoM_ref(CoM_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - CoM_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * q_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, data, t, dt, s), l_ref);
  cost->evalStageCostDerivatives(robot, data, t, dt, s, kkt_res);
  cost->evalStageCostHessian(robot, data, t, dt, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += dt * J_3d.transpose() * q_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_3d.transpose() * q_weight.asDiagonal() * J_3d;
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
  const Eigen::Vector3d q_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d CoM_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_weight(q_weight);
  cost->set_qf_weight(qf_weight);
  cost->set_qi_weight(qi_weight);
  cost->set_CoM_ref(CoM_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - CoM_ref;
  const double l_ref = 0.5 * q_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, t, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * qf_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * qf_weight.asDiagonal() * J_3d;
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
  const Eigen::Vector3d q_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d CoM_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_weight(q_weight);
  cost->set_qf_weight(qf_weight);
  cost->set_qi_weight(qi_weight);
  cost->set_CoM_ref(CoM_ref);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const Eigen::Vector3d q_diff = robot.CoM() - CoM_ref;
  const double l_ref = 0.5 * q_diff.transpose() * qi_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, data, t, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, data, t, s, kkt_res);
  cost->evalImpulseCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * qi_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * qi_weight.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


TEST_F(CoMCostTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(CoMCostTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}