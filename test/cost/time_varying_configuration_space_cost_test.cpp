#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/time_varying_configuration_space_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"

#include "idocp/utils/derivative_checker.hpp"


namespace idocp {

class TimeVaryingConfigurationSpaceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpulseCost(Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dt, t;
};


void TimeVaryingConfigurationSpaceCostTest::testStageCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qqq().setRandom();
  kkt_mat.Qvv().setRandom();
  kkt_mat.Qaa().setRandom();
  kkt_mat.Quu().setRandom();
  kkt_res.lq().setRandom();
  kkt_res.lv().setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu().setRandom();
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
  const double t_begin = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double t_end = t_begin + std::abs(Eigen::VectorXd::Random(1)[0]);
  const Eigen::VectorXd q0 = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v0 = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_q_weight(q_weight);
  cost->set_v_weight(v_weight);
  cost->set_a_weight(a_weight);
  cost->set_ref(robot, t_begin, t_end, q0, v0);
  Eigen::VectorXd q_ref = q0;
  Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(dimv);
  if (t > t_begin && t < t_end) {
    robot.integrateConfiguration(q0, v0, (t-t_begin), q_ref);
    v_ref = v0;
  }
  else if (t <= t_begin) {
    q_ref = q0;
    v_ref.setZero();
  }
  else {
    robot.integrateConfiguration(q0, v0, (t_end-t_begin), q_ref);
    v_ref.setZero();
  }
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
                            + (a_weight.array()*s.a.array()*s.a.array()).sum());
  EXPECT_DOUBLE_EQ(cost->computeStageCost(robot, data, t, dt, s), cost_ref);
  cost->computeStageCostDerivatives(robot, data, t, dt, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += dt * Jq_diff.transpose() * q_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += dt * q_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += dt * v_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.la += dt * a_weight.asDiagonal() * s.a;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->computeStageCostHessian(robot, data, t, dt, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += dt * Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += dt * q_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += dt * v_weight.asDiagonal();
  kkt_mat_ref.Qaa() += dt * a_weight.asDiagonal();
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


void TimeVaryingConfigurationSpaceCostTest::testTerminalCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qqq().setRandom();
  kkt_mat.Qvv().setRandom();
  kkt_res.lq().setRandom();
  kkt_res.lv().setRandom();
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  const double t_begin = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double t_end = t_begin + std::abs(Eigen::VectorXd::Random(1)[0]);
  const Eigen::VectorXd q0 = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v0 = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_qf_weight(qf_weight);
  cost->set_vf_weight(vf_weight);
  cost->set_ref(robot, t_begin, t_end, q0, v0);
  Eigen::VectorXd q_ref = q0;
  Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(dimv);
  if (t > t_begin && t < t_end) {
    robot.integrateConfiguration(q0, v0, (t-t_begin), q_ref);
    v_ref = v0;
  }
  else if (t <= t_begin) {
    q_ref = q0;
    v_ref.setZero();
  }
  else {
    robot.integrateConfiguration(q0, v0, (t_end-t_begin), q_ref);
    v_ref.setZero();
  }
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
  EXPECT_DOUBLE_EQ(cost->computeTerminalCost(robot, data, t, s), cost_ref);
  cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += qf_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += vf_weight.asDiagonal() * (s.v-v_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->computeTerminalCostHessian(robot, data, t, s, kkt_mat);
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


void TimeVaryingConfigurationSpaceCostTest::testImpulseCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  ImpulseSplitKKTMatrix kkt_mat(robot);
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_mat.Qqq().setRandom();
  kkt_mat.Qvv().setRandom();
  kkt_mat.Qdvdv().setRandom();
  kkt_res.lq().setRandom();
  kkt_res.lv().setRandom();
  kkt_res.ldv.setRandom();
  ImpulseSplitKKTMatrix kkt_mat_ref = kkt_mat;
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vi_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dvi_weight = Eigen::VectorXd::Random(dimv);
  const double t_begin = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double t_end = t_begin + std::abs(Eigen::VectorXd::Random(1)[0]);
  const Eigen::VectorXd q0 = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v0 = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost->useKinematics());
  cost->set_qi_weight(qi_weight);
  cost->set_vi_weight(vi_weight);
  cost->set_dvi_weight(dvi_weight);
  cost->set_ref(robot, t_begin, t_end, q0, v0);
  Eigen::VectorXd q_ref = q0;
  Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(dimv);
  if (t > t_begin && t < t_end) {
    robot.integrateConfiguration(q0, v0, (t-t_begin), q_ref);
    v_ref = v0;
  }
  else if (t <= t_begin) {
    q_ref = q0;
    v_ref.setZero();
  }
  else {
    robot.integrateConfiguration(q0, v0, (t_end-t_begin), q_ref);
    v_ref.setZero();
  }
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
  EXPECT_DOUBLE_EQ(cost->computeImpulseCost(robot, data, t, s), cost_ref);
  cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * qi_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += qi_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += vi_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.ldv += dvi_weight.asDiagonal() * s.dv;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->computeImpulseCostHessian(robot, data, t, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * qi_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += qi_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += vi_weight.asDiagonal();
  kkt_mat_ref.Qdvdv() += dvi_weight.asDiagonal();
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


TEST_F(TimeVaryingConfigurationSpaceCostTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(TimeVaryingConfigurationSpaceCostTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}