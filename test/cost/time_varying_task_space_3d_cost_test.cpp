#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"

namespace robotoc {
  
class TimeVaryingTaskSpace3DRef final : public TimeVaryingTaskSpace3DRefBase {
public:
  TimeVaryingTaskSpace3DRef(const Eigen::Vector3d& x3d0_ref, 
                            const Eigen::Vector3d& vx3d0_ref, 
                            const double t0, const double tf)
    : x3d0_ref_(x3d0_ref),
      vx3d0_ref_(vx3d0_ref),
      t0_(t0),
      tf_(tf) {
  }

  TimeVaryingTaskSpace3DRef() {}

  ~TimeVaryingTaskSpace3DRef() {}

  TimeVaryingTaskSpace3DRef(const TimeVaryingTaskSpace3DRef&) = default;

  TimeVaryingTaskSpace3DRef& operator=( 
      const TimeVaryingTaskSpace3DRef&) = default;

  TimeVaryingTaskSpace3DRef(
      TimeVaryingTaskSpace3DRef&&) noexcept = default;

  TimeVaryingTaskSpace3DRef& operator=(
      TimeVaryingTaskSpace3DRef&&) noexcept = default;

  void update_x3d_ref(const double t, Eigen::VectorXd& x3d_ref) const override {
    x3d_ref = x3d0_ref_ + (t-t0_) * vx3d0_ref_;
  }

  bool isActive(const double t) const override {
    if (t0_ <= t && t <= tf_)
      return true;
    else 
      return false;
  }

private:
  Eigen::Vector3d x3d0_ref_, vx3d0_ref_;
  double t0_, tf_;
};


class TimeVaryingTaskSpace3DCostTest : public ::testing::Test {
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

  void testStageCost(Robot& robot, const int frame_id) const;
  void testTerminalCost(Robot& robot, const int frame_id) const;
  void testImpulseCost(Robot& robot, const int frame_id) const;

  GridInfo grid_info, grid_info0, grid_infof;
  double dt, t, t0, tf;
};


void TimeVaryingTaskSpace3DCostTest::testStageCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3df_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3di_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx3d0_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TimeVaryingTaskSpace3DRef>(x3d0_ref, vx3d0_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, frame_id, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  cost->set_x3df_weight(x3df_weight);
  cost->set_x3di_weight(x3di_weight);
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

  const Eigen::Vector3d q_ref = x3d0_ref + (t-t0) * vx3d0_ref;
  const Eigen::Vector3d q_task = robot.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - q_ref;
  const double l_ref = dt * 0.5 * q_diff.transpose() * x3d_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, data, grid_info, s), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, data, grid_info, s, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot.frameRotation(frame_id) * J_6d.topRows(3);
  kkt_res_ref.lq() += dt * J_diff.transpose() * x3d_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dt * J_diff.transpose() * x3d_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


void TimeVaryingTaskSpace3DCostTest::testTerminalCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3df_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3di_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx3d0_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TimeVaryingTaskSpace3DRef>(x3d0_ref, vx3d0_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, frame_id, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  cost->set_x3df_weight(x3df_weight);
  cost->set_x3di_weight(x3di_weight);
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

  const Eigen::Vector3d q_ref = x3d0_ref + (t-t0) * vx3d0_ref;
  const Eigen::Vector3d q_task = robot.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - q_ref;
  const double l_ref = 0.5 * q_diff.transpose() * x3df_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, data, grid_info, s), l_ref);
  cost->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_res);
  cost->evalTerminalCostHessian(robot, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot.frameRotation(frame_id) * J_6d.topRows(3);
  kkt_res_ref.lq() += J_diff.transpose() * x3df_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_diff.transpose() * x3df_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void TimeVaryingTaskSpace3DCostTest::testImpulseCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::Vector3d x3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3df_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3di_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d x3d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx3d0_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TimeVaryingTaskSpace3DRef>(x3d0_ref, vx3d0_ref, t0, tf);
  auto cost = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, frame_id, ref);

  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_x3d_weight(x3d_weight);
  cost->set_x3df_weight(x3df_weight);
  cost->set_x3di_weight(x3di_weight);
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

  const Eigen::Vector3d q_ref = x3d0_ref + (t-t0) * vx3d0_ref;
  const Eigen::Vector3d q_task = robot.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - q_ref;
  const double l_ref = 0.5 * q_diff.transpose() * x3di_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost->evalImpulseCost(robot, impulse_status, data, grid_info, s), l_ref);
  cost->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, kkt_res);
  cost->evalImpulseCostHessian(robot, impulse_status, data, grid_info, s, kkt_mat);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot.frameRotation(frame_id) * J_6d.topRows(3);
  kkt_res_ref.lq() += J_diff.transpose() * x3di_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_diff.transpose() * x3di_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
}


TEST_F(TimeVaryingTaskSpace3DCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const int frame_id = robot.contactFrames()[0];
  testStageCost(robot, frame_id);
  testTerminalCost(robot, frame_id);
  testImpulseCost(robot, frame_id);
}


TEST_F(TimeVaryingTaskSpace3DCostTest, floatingBase) {
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