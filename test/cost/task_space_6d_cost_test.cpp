#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/cost/task_space_6d_cost.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"

#include "robotoc/utils/derivative_checker.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class TestTaskSpace6DRef final : public TaskSpace6DRefBase {
public:
  TestTaskSpace6DRef(const Eigen::Vector3d& x6d0_ref, 
                            const Eigen::Vector3d& vx6d_ref, 
                            const double t0, const double tf)
    : x6d0_ref_(x6d0_ref),
      vx6d_ref_(vx6d_ref),
      t0_(t0),
      tf_(tf),
      rotm_(Eigen::Matrix3d::Identity()) {
  }

  TestTaskSpace6DRef() {}

  ~TestTaskSpace6DRef() {}

  TestTaskSpace6DRef(const TestTaskSpace6DRef&) = default;

  TestTaskSpace6DRef& operator=(const TestTaskSpace6DRef&) = default;

  TestTaskSpace6DRef(TestTaskSpace6DRef&&) noexcept = default;

  TestTaskSpace6DRef& operator=(TestTaskSpace6DRef&&) noexcept = default;

  void updateRef(const GridInfo& grid_info, SE3& x6d_ref) const override {
    x6d_ref = SE3(rotm_, (x6d0_ref_+(grid_info.t-t0_)*vx6d_ref_));
  }

  bool isActive(const GridInfo& grid_info) const override {
    if (t0_ <= grid_info.t && grid_info.t <= tf_)
      return true;
    else 
      return false;
  }


private:
  Eigen::Vector3d x6d0_ref_, vx6d_ref_;
  Eigen::Matrix3d rotm_;
  double t0_, tf_;
};



class TaskSpace6DCostTest : public ::testing::Test {
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

  void testStageCostConstRef(Robot& robot, const int frame_id) const;
  void testTerminalCostConstRef(Robot& robot, const int frame_id) const;
  void testImpactCostConstRef(Robot& robot, const int frame_id) const;
  void testStageCost(Robot& robot, const int frame_id) const;
  void testTerminalCost(Robot& robot, const int frame_id) const;
  void testImpactCost(Robot& robot, const int frame_id) const;

  GridInfo grid_info, grid_info0, grid_infof;
  double dt, t, t0, tf;
};


TEST_F(TaskSpace6DCostTest, defaultConstructor) {
  EXPECT_NO_THROW(
    auto cost = std::make_shared<TaskSpace6DCost>();
  );
}


void TaskSpace6DCostTest::testStageCostConstRef(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, robot.frameName(frame_id));
  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  cost->set_const_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight.asDiagonal() * diff_6d;
  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_info, s, data), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
}


void TaskSpace6DCostTest::testTerminalCostConstRef(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  cost->set_const_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight_terminal.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_info, s, data), l_ref);
  cost->evalTerminalCostDerivatives(robot, grid_info, s, data, kkt_res);
  cost->evalTerminalCostHessian(robot, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight_terminal.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight_terminal.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void TaskSpace6DCostTest::testImpactCostConstRef(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const SE3 x6d_ref = SE3::Random();
  const Eigen::Vector3d position_ref = x6d_ref.translation();
  const Eigen::Matrix3d rotation_ref = x6d_ref.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  cost->set_const_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight_impact.asDiagonal() * diff_6d;
  const auto impact_status = robot.createImpactStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_info, s, data), l_ref);
  cost->evalImpactCostDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  cost->evalImpactCostHessian(robot, impact_status, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight_impact.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight_impact.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpactCostDerivatives(cost));
}


void TaskSpace6DCostTest::testStageCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::Vector3d x6d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx6d_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TestTaskSpace6DRef>(x6d0_ref, vx6d_ref, t0, tf);
  auto cost = std::make_shared<TaskSpace6DCost>(robot, robot.frameName(frame_id), ref);

  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

  const auto contact_status = robot.createContactStatus();
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_info0, s, data), 0);
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_infof, s, data), 0);
  cost->evalStageCostDerivatives(robot, contact_status, grid_info0, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostDerivatives(robot, contact_status, grid_infof, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalStageCostHessian(robot, contact_status, grid_info0, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalStageCostHessian(robot, contact_status, grid_infof, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d q_ref = x6d0_ref + (t-t0) * vx6d_ref;
  const SE3 x6d_ref = SE3(Eigen::Matrix3d::Identity(), q_ref);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->evalStageCost(robot, contact_status, grid_info, s, data), l_ref);
  cost->evalStageCostDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  cost->evalStageCostHessian(robot, contact_status, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight.asDiagonal() * J66_J_6d;
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
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::Vector3d x6d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx6d_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TestTaskSpace6DRef>(x6d0_ref, vx6d_ref, t0, tf);
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id, ref);

  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);

  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_info0, s, data), 0);
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_infof, s, data), 0);
  cost->evalTerminalCostDerivatives(robot, grid_info0, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostDerivatives(robot, grid_infof, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalTerminalCostHessian(robot, grid_info0, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalTerminalCostHessian(robot, grid_infof, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d q_ref = x6d0_ref + (t-t0) * vx6d_ref;
  const SE3 x6d_ref = SE3(Eigen::Matrix3d::Identity(), q_ref);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight_terminal.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->evalTerminalCost(robot, grid_info, s, data), l_ref);
  cost->evalTerminalCostDerivatives(robot, grid_info, s, data, kkt_res);
  cost->evalTerminalCostHessian(robot, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight_terminal.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight_terminal.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
}


void TaskSpace6DCostTest::testImpactCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_terminal = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd weight_impact = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::Vector3d x6d0_ref = Eigen::Vector3d::Random();
  const Eigen::Vector3d vx6d_ref = Eigen::Vector3d::Random();
  auto ref = std::make_shared<TestTaskSpace6DRef>(x6d0_ref, vx6d_ref, t0, tf);
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id, ref);

  CostFunctionData data(robot);
  cost->set_weight(weight.tail(3), weight.head(3));
  cost->set_weight_terminal(weight_terminal.tail(3), weight_terminal.head(3));
  cost->set_weight_impact(weight_impact.tail(3), weight_impact.head(3));
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);

  const auto impact_status = robot.createImpactStatus();
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_info0, s, data), 0);
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_infof, s, data), 0);
  cost->evalImpactCostDerivatives(robot, impact_status, grid_info0, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostDerivatives(robot, impact_status, grid_infof, s, data, kkt_res);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost->evalImpactCostHessian(robot, impact_status, grid_info0, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost->evalImpactCostHessian(robot, impact_status, grid_infof, s, data, kkt_mat);
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));

  const Eigen::Vector3d q_ref = x6d0_ref + (t-t0) * vx6d_ref;
  const SE3 x6d_ref = SE3(Eigen::Matrix3d::Identity(), q_ref);
  const SE3 placement = robot.framePlacement(frame_id);
  const SE3 diff_x6d = x6d_ref.inverse() * placement;
  const Eigen::VectorXd diff_6d = Log6Map(diff_x6d);
  const double l_ref = 0.5 * diff_6d.transpose() * weight_impact.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->evalImpactCost(robot, impact_status, grid_info, s, data), l_ref);
  cost->evalImpactCostDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  cost->evalImpactCostHessian(robot, impact_status, grid_info, s, data, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  computeJLog6Map(diff_x6d, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * weight_impact.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * weight_impact.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-02);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpactCostDerivatives(cost));
}


TEST_F(TaskSpace6DCostTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const int frame_id = robot.contactFrames()[0];
  testStageCostConstRef(robot, frame_id);
  testTerminalCostConstRef(robot, frame_id);
  testImpactCostConstRef(robot, frame_id);
  testStageCost(robot, frame_id);
  testTerminalCost(robot, frame_id);
  testImpactCost(robot, frame_id);
}


TEST_F(TaskSpace6DCostTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const std::vector<int> frames = robot.contactFrames();
  for (const auto frame_id : frames) {
    testStageCostConstRef(robot, frame_id);
    testTerminalCostConstRef(robot, frame_id);
    testImpactCostConstRef(robot, frame_id);
    testStageCost(robot, frame_id);
    testTerminalCost(robot, frame_id);
    testImpactCost(robot, frame_id);
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}