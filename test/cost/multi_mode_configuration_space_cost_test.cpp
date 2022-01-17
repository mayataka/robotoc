#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/multi_mode_configuration_space_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class MultiModeConfigurationSpaceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    t = grid_info.t;
    dt = grid_info.dt;
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpulseCost(Robot& robot) const;

  GridInfo grid_info;
  double dt, t;
};


void MultiModeConfigurationSpaceCostTest::testStageCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qaa.setRandom();
  kkt_mat.Quu.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimu);
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimu);
  auto cost = std::make_shared<MultiModeConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_q_weight(q_weight, 0);
  cost->set_v_weight(v_weight, 0);
  cost->set_a_weight(a_weight, 0);
  cost->set_u_weight(u_weight, 0);
  cost->set_q_ref(q_ref, 0);
  cost->set_q_ref(q_ref, {0, 1, 2, 3});
  cost->set_v_ref(v_ref, 0);
  cost->set_u_ref(u_ref, 0);
  const SplitSolution s = SplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double cost_ref = 0.5 * dt 
                           * ((q_weight.array()*q_diff.array()*q_diff.array()).sum()
                            + (v_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum()
                            + (a_weight.array()*s.a.array()*s.a.array()).sum()
                            + (u_weight.array()* (s.u-u_ref).array()*(s.u-u_ref).array()).sum());
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_info, s), cost_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += dt * Jq_diff.transpose() * q_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += dt * q_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += dt * v_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.la += dt * a_weight.asDiagonal() * s.a;
  kkt_res_ref.lu += dt * u_weight.asDiagonal() * (s.u-u_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostHessian(robot, contact_status, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += dt * Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += dt * q_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += dt * v_weight.asDiagonal();
  kkt_mat_ref.Qaa += dt * a_weight.asDiagonal();
  kkt_mat_ref.Quu += dt * u_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
  if (robot.hasFloatingBase()) {
    // This is due to Gauss-Newton Hessian approximation.
    EXPECT_FALSE(derivative_checker.checkSecondOrderStageCostDerivatives(cost));
  }
  else {
    EXPECT_TRUE(derivative_checker.checkSecondOrderStageCostDerivatives(cost));
  }
}


void MultiModeConfigurationSpaceCostTest::testTerminalCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_res.lx.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<MultiModeConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_qf_weight(qf_weight);
  cost->set_vf_weight(vf_weight);
  cost->set_q_ref(q_ref, 0);
  cost->set_v_ref(v_ref, 0);
  const SplitSolution s = SplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double cost_ref = 0.5 * ((qf_weight.array()* q_diff.array()*q_diff.array()).sum()
                                + (vf_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), cost_ref);
  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += qf_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += vf_weight.asDiagonal() * (s.v-v_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * qf_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += qf_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += vf_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
  if (robot.hasFloatingBase()) {
    // This is due to Gauss-Newton Hessian approximation.
    EXPECT_FALSE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
  }
  else {
    EXPECT_TRUE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
  }
}


void MultiModeConfigurationSpaceCostTest::testImpulseCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  ImpulseSplitKKTMatrix kkt_mat(robot);
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qdvdv.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.ldv.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vi_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dvi_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<MultiModeConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_qi_weight(qi_weight);
  cost->set_vi_weight(vi_weight);
  cost->set_dvi_weight(dvi_weight);
  cost->set_q_ref(q_ref, 0);
  cost->set_v_ref(v_ref, 0);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double cost_ref = 0.5 * ((qi_weight.array()* q_diff.array()*q_diff.array()).sum()
                                + (vi_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum()
                                + (dvi_weight.array()* s.dv.array()*s.dv.array()).sum());
  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_info, s), cost_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * qi_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += qi_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += vi_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.ldv += dvi_weight.asDiagonal() * s.dv;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * qi_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += qi_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += vi_weight.asDiagonal();
  kkt_mat_ref.Qdvdv += dvi_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
  if (robot.hasFloatingBase()) {
    // This is due to Gauss-Newton Hessian approximation.
    EXPECT_FALSE(derivative_checker.checkSecondOrderImpulseCostDerivatives(cost));
  }
  else {
    EXPECT_TRUE(derivative_checker.checkSecondOrderImpulseCostDerivatives(cost));
  }
}


TEST_F(MultiModeConfigurationSpaceCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(MultiModeConfigurationSpaceCostTest, floatingBase) {
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