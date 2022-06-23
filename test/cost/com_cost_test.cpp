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

class TestCoMRef final : public CoMRefBase {
public:
  TestCoMRef(const Eigen::Vector3d& com0_ref, const Eigen::Vector3d& vcom_ref, 
             const double t0, const double tf)
    : com0_ref_(com0_ref),
      vcom_ref_(vcom_ref),
      t0_(t0),
      tf_(tf) {
  }

  TestCoMRef() {}

  ~TestCoMRef() {}

  TestCoMRef(const TestCoMRef&) = default;

  TestCoMRef& operator=(const TestCoMRef&) = default;

  TestCoMRef(TestCoMRef&&) noexcept = default;

  TestCoMRef& operator=(TestCoMRef&&) noexcept = default;

  void updateRef(const GridInfo& grid_info, 
                 Eigen::VectorXd& ref) const override {
    ref = com0_ref_ + (grid_info.t-t0_) * vcom_ref_;
  }

  bool isActive(const GridInfo& grid_info) const override {
    if (t0_ <= grid_info.t && grid_info.t <= tf_)
      return true;
    else 
      return false;
  }

private:
  Eigen::Vector3d com0_ref_, vcom_ref_;
  double t0_, tf_;
};


class CoMCostTest : public ::testing::Test {
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
  void testImpulseCostConstRef(Robot& robot) const;
  void testStageCost(Robot& robot) const;
  void testTerminalCost(Robot& robot) const;
  void testImpulseCost(Robot& robot) const;

  GridInfo grid_info, grid_info0, grid_infof;
  double dt, t, t0, tf;
};


void CoMCostTest::testStageCostConstRef(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d const_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  cost->set_const_ref(const_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - const_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * weight.asDiagonal() * q_diff;
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_info, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_info, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += dt * J_3d.transpose() * weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_3d.transpose() * weight.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


void CoMCostTest::testTerminalCostConstRef(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d const_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  cost->set_const_ref(const_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_diff = robot.CoM() - const_ref;
  const double l_ref = 0.5 * q_diff.transpose() * weight_terminal.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * weight_terminal.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * weight_terminal.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void CoMCostTest::testImpulseCostConstRef(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d const_ref = Eigen::Vector3d::Random();
  auto cost = std::make_shared<CoMCost>(robot);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  cost->set_const_ref(const_ref);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const Eigen::Vector3d q_diff = robot.CoM() - const_ref;
  const double l_ref = 0.5 * q_diff.transpose() * weight_impulse.asDiagonal() * q_diff;
  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_info, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * weight_impulse.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * weight_impulse.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


void CoMCostTest::testStageCost(Robot& robot) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TestCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<CoMCost>(robot, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

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

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_info, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_info, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += dt * J_3d.transpose() * weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_3d.transpose() * weight.asDiagonal() * J_3d;
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
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TestCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<CoMCost>(robot, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

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

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = 0.5 * q_diff.transpose() * weight_terminal.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * weight_terminal.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * weight_terminal.asDiagonal() * J_3d;
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
  const Eigen::Vector3d weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_terminal = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d weight_impulse = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TestCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<CoMCost>(robot);
  cost->set_ref(ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_weight(weight);
  cost->set_weight_terminal(weight_terminal);
  cost->set_weight_impulse(weight_impulse);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);

  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_info0, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_infof, s), 0);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info0, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_infof, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_info0, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_infof, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = 0.5 * q_diff.transpose() * weight_impulse.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_info, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_3d = Eigen::MatrixXd::Zero(3, dimv);
  robot.getCoMJacobian(J_3d);
  kkt_res_ref.lq() += J_3d.transpose() * weight_impulse.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_3d.transpose() * weight_impulse.asDiagonal() * J_3d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


TEST_F(CoMCostTest, defaultConstructor) {
  EXPECT_NO_THROW(
    auto cost = std::make_shared<CoMCost>();
  );
}


TEST_F(CoMCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCostConstRef(robot);
  testTerminalCostConstRef(robot);
  testImpulseCostConstRef(robot);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}


TEST_F(CoMCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  testStageCostConstRef(robot);
  testTerminalCostConstRef(robot);
  testImpulseCostConstRef(robot);
  testStageCost(robot);
  testTerminalCost(robot);
  testImpulseCost(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}