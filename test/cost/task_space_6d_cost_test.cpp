#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/cost/task_space_6d_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class TaskSpace6DCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot, const int frame_id) const;
  void testTerminalCost(Robot& robot, const int frame_id) const;
  void testImpulseCost(Robot& robot, const int frame_id) const;

  double dt, t;
};


void TaskSpace6DCostTest::testStageCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd x6d_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6df_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6di_weight = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x6d_weight(x6d_weight.tail(3), x6d_weight.head(3));
  cost->set_x6df_weight(x6df_weight.tail(3), x6df_weight.head(3));
  cost->set_x6di_weight(x6di_weight.tail(3), x6di_weight.head(3));
  cost->set_x6d_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = dt * 0.5 * diff_6d.transpose() * x6d_weight.asDiagonal() * diff_6d;
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, t, dt, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, t, dt, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += dt * J66_J_6d.transpose() * x6d_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += dt * J66_J_6d.transpose() * x6d_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


void TaskSpace6DCostTest::testTerminalCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd x6d_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6df_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6di_weight = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x6d_weight(x6d_weight.tail(3), x6d_weight.head(3));
  cost->set_x6df_weight(x6df_weight.tail(3), x6df_weight.head(3));
  cost->set_x6di_weight(x6di_weight.tail(3), x6di_weight.head(3));
  cost->set_x6d_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * x6df_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, t, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * x6df_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * x6df_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void TaskSpace6DCostTest::testImpulseCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd x6d_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6df_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd x6di_weight = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x6d_weight(x6d_weight.tail(3), x6d_weight.head(3));
  cost->set_x6df_weight(x6df_weight.tail(3), x6df_weight.head(3));
  cost->set_x6di_weight(x6di_weight.tail(3), x6di_weight.head(3));
  cost->set_x6d_ref(position_ref, rotation_ref);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * x6di_weight.asDiagonal() * diff_6d;
  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, t, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, t, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, t, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * x6di_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * x6di_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


TEST_F(TaskSpace6DCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const int frame_id = robot.contactFrames()[0];
  testStageCost(robot, frame_id);
  testTerminalCost(robot, frame_id);
  testImpulseCost(robot, frame_id);
}


TEST_F(TaskSpace6DCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const std::vector<int> frames = robot.contactFrames();
  for (const auto frame_id : frames) {
    testStageCost(robot, frame_id);
    testTerminalCost(robot, frame_id);
    testImpulseCost(robot, frame_id);
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}