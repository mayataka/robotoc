#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_com_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class TimeVaryingCoMRef final : public TimeVaryingCoMRefBase {
public:
  TimeVaryingCoMRef(const Eigen::Vector3d& com0_ref, const Eigen::Vector3d& vcom_ref, 
                    const double t0, const double tf)
    : com0_ref_(com0_ref),
      vcom_ref_(vcom_ref),
      t0_(t0),
      tf_(tf) {
  }

  TimeVaryingCoMRef() {}

  ~TimeVaryingCoMRef() {}

  TimeVaryingCoMRef(const TimeVaryingCoMRef&) = default;

  TimeVaryingCoMRef& operator=(const TimeVaryingCoMRef&) = default;

  TimeVaryingCoMRef(TimeVaryingCoMRef&&) noexcept = default;

  TimeVaryingCoMRef& operator=(TimeVaryingCoMRef&&) noexcept = default;

  void update_com_ref(const double t, Eigen::VectorXd& com_ref) const override {
    com_ref = com0_ref_ + (t-t0_) * vcom_ref_;
  }

  bool isActive(const double t) const override {
    if (t0_ <= t && t <= tf_)
      return true;
    else 
      return false;
  }

private:
  Eigen::Vector3d com0_ref_, vcom_ref_;
  double t0_, tf_;
};


class TimeVaryingCoMCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    t0 = t - std::abs(Eigen::VectorXd::Random(1)[0]);
    tf = t + std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot, const int frame_id) const;
  void testTerminalCost(Robot& robot, const int frame_id) const;
  void testImpulseCost(Robot& robot, const int frame_id) const;

  double dt, t, t0, tf;
};


void TimeVaryingCoMCostTest::testStageCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TimeVaryingCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingCoMCost>(robot, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, t0-dt, dt, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, tf+dt, dt, s), 0);
  cost->evalStageCostDerivatives(robot, contact_status, data, t0-dt, dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostDerivatives(robot, contact_status, data, tf+dt, dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostHessian(robot, contact_status, data, t0-dt, dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalStageCostHessian(robot, contact_status, data, tf+dt, dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * com_weight.asDiagonal() * q_diff;
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


void TimeVaryingCoMCostTest::testTerminalCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TimeVaryingCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingCoMCost>(robot, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, t0-dt, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, tf+dt, s), 0);
  cost->evalTerminalCostDerivatives(robot, data, t0-dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostDerivatives(robot, data, tf+dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, data, t0-dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalTerminalCostHessian(robot, data, tf+dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
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


void TimeVaryingCoMCostTest::testImpulseCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d com_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d comi_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d com0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vcom_ref = Eigen::Vector3d::Random();

  auto ref = std::make_shared<TimeVaryingCoMRef>(com0_ref, vcom_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingCoMCost>(robot, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_com_weight(com_weight);
  cost->set_comf_weight(comf_weight);
  cost->set_comi_weight(comi_weight);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);

  const auto impulse_status = robot.createImpulseStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, t0-dt, s), 0);
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, tf+dt, s), 0);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, t0-dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, tf+dt, s, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpulseCostHessian(robot, impulse_status, data, t0-dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalImpulseCostHessian(robot, impulse_status, data, tf+dt, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d com_ref = com0_ref + (t-t0) * vcom_ref;
  const Eigen::Vector3d q_diff = robot.CoM() - com_ref;
  const double l_ref = 0.5 * q_diff.transpose() * comi_weight.asDiagonal() * q_diff;
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


TEST_F(TimeVaryingCoMCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const int frame_id = robot.contactFrames()[0];
  testStageCost(robot, frame_id);
  testTerminalCost(robot, frame_id);
  testImpulseCost(robot, frame_id);
}


TEST_F(TimeVaryingCoMCostTest, floatingBase) {
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