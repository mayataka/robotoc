#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/swing_foot_cost.hpp"
#include "robotoc/cost/trotting_swing_foot_ref.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class SwingFootCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot, const int contact_index) const;
  void testTerminalCost(Robot& robot, const int contact_index) const;
  void testImpulseCost(Robot& robot, const int contact_index) const;

  double dt, t;
};


TEST_F(SwingFootCostTest, testStageCost) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const int contact_index = 0;
  const int x_ref_foot_contact_index = 1;
  const int y_ref_foot_contact_index = 2;
  const double step_length = 1.5;
  const double step_height = 1.0;
  auto ref = std::make_shared<TrottingSwingFootRef>(contact_index,
                                                    x_ref_foot_contact_index, 
                                                    y_ref_foot_contact_index,
                                                    step_length, step_height);
  auto cost = std::make_shared<SwingFootCost>(robot, contact_index, ref);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  cost->set_x3d_ref(ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random());
  }
  Eigen::Vector3d x3d_ref;
  const double xdiff = contact_status.contactPosition(contact_index)[0] 
                        - contact_status.contactPosition(x_ref_foot_contact_index)[0];
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (std::abs(xdiff) < eps) {
    x3d_ref[0] = contact_status.contactPosition(x_ref_foot_contact_index)[0];
    x3d_ref[0] += 0.25 * step_length;
  }
  else {
    x3d_ref[0] = contact_status.contactPosition(x_ref_foot_contact_index)[0];
  }
  x3d_ref[1] = contact_status.contactPosition(y_ref_foot_contact_index)[1];
  x3d_ref[2] = step_height;
  const int frame_id = robot.contactFrames()[contact_index];
  const Eigen::Vector3d q_task = robot.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - x3d_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * x3d_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, t, dt, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, t, dt, s, kkt_mat);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot.frameRotation(frame_id) * J_6d.topRows(3);
  kkt_res_ref.lq() += dt * J_diff.transpose() * x3d_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_diff.transpose() * x3d_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    contact_status.activateContact(i);
  }
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, t, dt, s), 0.);
  cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, t, dt, s, kkt_mat);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


TEST_F(SwingFootCostTest, testTerminalCost) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const int contact_index = 0;
  const int x_ref_foot_contact_index = 1;
  const int y_ref_foot_contact_index = 2;
  const double step_length = 1.5;
  const double step_height = 1.0;
  auto ref = std::make_shared<TrottingSwingFootRef>(contact_index,
                                                    x_ref_foot_contact_index, 
                                                    y_ref_foot_contact_index,
                                                    step_length, step_height);
  auto cost = std::make_shared<SwingFootCost>(robot, contact_index, ref);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const double l_ref = 0;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, t, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


TEST_F(SwingFootCostTest, testImpulseCost) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const int contact_index = 0;
  const int x_ref_foot_contact_index = 1;
  const int y_ref_foot_contact_index = 2;
  const double step_length = 1.5;
  const double step_height = 1.0;
  auto ref = std::make_shared<TrottingSwingFootRef>(contact_index,
                                                    x_ref_foot_contact_index, 
                                                    y_ref_foot_contact_index,
                                                    step_length, step_height);
  auto cost = std::make_shared<SwingFootCost>(robot, contact_index, ref);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const double l_ref = 0.0;
  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, t, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, t, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}