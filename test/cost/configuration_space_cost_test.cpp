#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class TestConfigurationSpaceRef final : public ConfigurationSpaceRefBase {
public:
  TestConfigurationSpaceRef(const Eigen::VectorXd& q0_ref, 
                              const Eigen::VectorXd& v_ref,
                              const double t0, const double tf)
    : q0_ref_(q0_ref),
      v_ref_(v_ref),
      t0_(t0),
      tf_(tf) {
  }

  TestConfigurationSpaceRef() {}

  ~TestConfigurationSpaceRef() {}

  TestConfigurationSpaceRef(const TestConfigurationSpaceRef&) = default;

  TestConfigurationSpaceRef& operator=( 
      const TestConfigurationSpaceRef&) = default;

  TestConfigurationSpaceRef(
      TestConfigurationSpaceRef&&) noexcept = default;

  TestConfigurationSpaceRef& operator=(
      TestConfigurationSpaceRef&&) noexcept = default;

  void updateRef(const Robot& robot, const GridInfo& grid_info,
                 Eigen::VectorXd& q_ref) const override {
    robot.integrateConfiguration(q0_ref_, v_ref_, (grid_info.t-t0_), q_ref);
  }

  bool isActive(const GridInfo& grid_info) const override {
    if (t0_ <= grid_info.t && grid_info.t <= tf_)
      return true;
    else 
      return false;
  }

private:
  Eigen::VectorXd q0_ref_, v_ref_;
  double t0_, tf_;
};

class ConfigurationSpaceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    grid_info0 = grid_info;
    grid_infof = grid_info;
    dt = grid_info0.dt;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    t0 = t - std::abs(Eigen::VectorXd::Random(1)[0]);
    tf = t + std::abs(Eigen::VectorXd::Random(1)[0]);
    grid_info.t = t;
    grid_info0.t = t0 - dt;
    grid_infof.t = tf + dt;
  }

  virtual void TearDown() {
  }

  void testStageCostConstRef(Robot& robot) const;
  void testTerminalCostConstRef(Robot& robot) const;
  void testImpactCostConstRef(Robot& robot) const;
  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpactCost(Robot& robot) const;

  GridInfo grid_info, grid_info0, grid_infof;
  double dt, t, t0, tf;
};


void ConfigurationSpaceCostTest::testStageCostConstRef(Robot& robot) const {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv).cwiseAbs(); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimu).cwiseAbs();
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimu);
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  cost->set_q_weight(q_weight);
  cost->set_v_weight(v_weight);
  cost->set_a_weight(a_weight);
  cost->set_u_weight(u_weight);
  cost->set_q_ref(q_ref);
  cost->set_v_ref(v_ref);
  cost->set_u_ref(u_ref);
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


void ConfigurationSpaceCostTest::testTerminalCostConstRef(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_res.lx.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight_terminal = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd v_weight_terminal = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  cost->set_q_weight_terminal(q_weight_terminal);
  cost->set_v_weight_terminal(v_weight_terminal);
  cost->set_q_ref(q_ref);
  cost->set_v_ref(v_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double cost_ref = 0.5 * ((q_weight_terminal.array()* q_diff.array()*q_diff.array()).sum()
                                + (v_weight_terminal.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), cost_ref);
  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * q_weight_terminal.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += q_weight_terminal.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += v_weight_terminal.asDiagonal() * (s.v-v_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * q_weight_terminal.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += q_weight_terminal.asDiagonal();
  }
  kkt_mat_ref.Qvv() += v_weight_terminal.asDiagonal();
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


void ConfigurationSpaceCostTest::testImpactCostConstRef(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qdvdv.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.ldv.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight_impact = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd v_weight_impact = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd dv_weight_impact = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot);
  CostFunctionData data(robot);
  cost->set_q_weight_impact(q_weight_impact);
  cost->set_v_weight_impact(v_weight_impact);
  cost->set_dv_weight_impact(dv_weight_impact);
  cost->set_q_ref(q_ref);
  cost->set_v_ref(v_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double cost_ref = 0.5 * ((q_weight_impact.array()* q_diff.array()*q_diff.array()).sum()
                                + (v_weight_impact.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum()
                                + (dv_weight_impact.array()* s.dv.array()*s.dv.array()).sum());
  const auto impact_status = robot.createImpactStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, data, grid_info, s), cost_ref);
  cost->evalImpactCostDerivatives(robot, impact_status, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * q_weight_impact.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += q_weight_impact.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += v_weight_impact.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.ldv += dv_weight_impact.asDiagonal() * s.dv;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostHessian(robot, impact_status, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * q_weight_impact.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += q_weight_impact.asDiagonal();
  }
  kkt_mat_ref.Qvv() += v_weight_impact.asDiagonal();
  kkt_mat_ref.Qdvdv += dv_weight_impact.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpactCostDerivatives(cost));
  if (robot.hasFloatingBase()) {
    // This is due to Gauss-Newton Hessian approximation.
    EXPECT_FALSE(derivative_checker.checkSecondOrderImpactCostDerivatives(cost));
  }
  else {
    EXPECT_TRUE(derivative_checker.checkSecondOrderImpactCostDerivatives(cost));
  }
}


void ConfigurationSpaceCostTest::testStageCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd q0_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto ref = std::make_shared<TestConfigurationSpaceRef>(q0_ref, v_ref, t0, tf);
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot, ref);
  CostFunctionData data(robot);
  cost->set_q_weight(q_weight);

  const auto s = SplitSolution::Random(robot);
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_info0, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_infof, s), 0);
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_info0, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_infof, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostHessian(robot, contact_status, data, grid_info0, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalStageCostHessian(robot, contact_status, data, grid_infof, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot.integrateConfiguration(q0_ref, v_ref, (t-t0), q_ref);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  const double cost_ref = 0.5 * dt * (q_weight.array()*q_diff.array()*q_diff.array()).sum();
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
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostHessian(robot, contact_status, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += dt * Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += dt * q_weight.asDiagonal();
  }
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
  if (!robot.hasFloatingBase()) {
    EXPECT_TRUE(derivative_checker.checkSecondOrderStageCostDerivatives(cost));
  }
}


void ConfigurationSpaceCostTest::testTerminalCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight_terminal = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd q0_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto ref = std::make_shared<TestConfigurationSpaceRef>(q0_ref, v_ref, t0, tf);
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot, ref);
  CostFunctionData data(robot);
  cost->set_q_weight_terminal(q_weight_terminal);

  const auto s = SplitSolution::Random(robot);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info0, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_infof, s), 0);
  cost->evalTerminalCostDerivatives(robot, data, grid_info0, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostDerivatives(robot, data, grid_infof, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, data, grid_info0, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalTerminalCostHessian(robot, data, grid_infof, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot.integrateConfiguration(q0_ref, v_ref, (t-t0), q_ref);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  const double cost_ref = 0.5 * (q_weight_terminal.array()*q_diff.array()*q_diff.array()).sum();
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), cost_ref);

  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * q_weight_terminal.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += q_weight_terminal.asDiagonal() * (s.q-q_ref);
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * q_weight_terminal.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += q_weight_terminal.asDiagonal();
  }
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
  if (!robot.hasFloatingBase()) {
    EXPECT_TRUE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
  }
}


void ConfigurationSpaceCostTest::testImpactCost(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight_impact = Eigen::VectorXd::Random(dimv).cwiseAbs();
  const Eigen::VectorXd q0_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  auto ref = std::make_shared<TestConfigurationSpaceRef>(q0_ref, v_ref, t0, tf);
  auto cost = std::make_shared<ConfigurationSpaceCost>(robot, ref);
  CostFunctionData data(robot);
  cost->set_q_weight_impact(q_weight_impact);

  const auto s = SplitSolution::Random(robot);
  const auto impact_status = robot.createImpactStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, data, grid_info0, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, data, grid_infof, s), 0);
  cost->evalImpactCostDerivatives(robot, impact_status, data, grid_info0, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostDerivatives(robot, impact_status, data, grid_infof, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostHessian(robot, impact_status, data, grid_info0, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalImpactCostHessian(robot, impact_status, data, grid_infof, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot.integrateConfiguration(q0_ref, v_ref, (t-t0), q_ref);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  robot.subtractConfiguration(s.q, q_ref, q_diff);
  const double cost_ref = 0.5 * (q_weight_impact.array()*q_diff.array()*q_diff.array()).sum();
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, data, grid_info, s), cost_ref);

  cost->evalImpactCostDerivatives(robot, impact_status, data, grid_info, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * q_weight_impact.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += q_weight_impact.asDiagonal() * (s.q-q_ref);
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostHessian(robot, impact_status, data, grid_info, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * q_weight_impact.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += q_weight_impact.asDiagonal();
  }
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpactCostDerivatives(cost));
  if (!robot.hasFloatingBase()) {
    EXPECT_TRUE(derivative_checker.checkSecondOrderImpactCostDerivatives(cost));
  }
}


TEST_F(ConfigurationSpaceCostTest, defaultConstructor) {
  EXPECT_NO_THROW(
    auto cost = std::make_shared<ConfigurationSpaceCost>();
  );
}


TEST_F(ConfigurationSpaceCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCostConstRef(robot);
  testTerminalCostConstRef(robot);
  testImpactCostConstRef(robot);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpactCost(robot);
}


TEST_F(ConfigurationSpaceCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  testStageCostConstRef(robot);
  testTerminalCostConstRef(robot);
  testImpactCostConstRef(robot);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpactCost(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}