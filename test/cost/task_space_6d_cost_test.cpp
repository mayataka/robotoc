#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/task_space_6d_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"

#include "idocp/utils/derivative_checker.hpp"


namespace idocp {

class TaskSpace6DCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot, const int frame_id) const;
  void testTerminalCost(Robot& robot, const int frame_id) const;
  void testImpulseCost(Robot& robot, const int frame_id) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau, t;
};


void TaskSpace6DCostTest::testStageCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(6).array().abs();
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost->set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost->set_qi_6d_weight(qi_weight.tail(3), qi_weight.head(3));
  cost->set_q_6d_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();
  const double l_ref = dtau * 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->computeStageCost(robot, data, t, dtau, s), l_ref);
  cost->computeStageCostDerivatives(robot, data, t, dtau, s, kkt_res);
  cost->computeStageCostHessian(robot, data, t, dtau, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::Jlog6(diff_SE3, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-03);
  EXPECT_TRUE(derivative_checker.checkFirstOrderStageCostDerivatives(cost));
  // This is due to Gauss-Newton Hessian approximation.
  EXPECT_FALSE(derivative_checker.checkSecondOrderStageCostDerivatives(cost));
}


void TaskSpace6DCostTest::testTerminalCost(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(6).array().abs();
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost->set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost->set_qi_6d_weight(qi_weight.tail(3), qi_weight.head(3));
  cost->set_q_6d_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();
  const double l_ref = 0.5 * diff_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->computeTerminalCost(robot, data, t, s), l_ref);
  cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeTerminalCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::Jlog6(diff_SE3, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-03);
  EXPECT_TRUE(derivative_checker.checkFirstOrderTerminalCostDerivatives(cost));
  // This is due to Gauss-Newton Hessian approximation.
  EXPECT_FALSE(derivative_checker.checkSecondOrderTerminalCostDerivatives(cost));
}


void TaskSpace6DCostTest::testImpulseCost(Robot& robot, const int frame_id) const {
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(6).array().abs();
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(6).array().abs();
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  auto cost = std::make_shared<TaskSpace6DCost>(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost->useKinematics());
  cost->set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost->set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost->set_qi_6d_weight(qi_weight.tail(3), qi_weight.head(3));
  cost->set_q_6d_ref(position_ref, rotation_ref);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v);
  const pinocchio::SE3 placement = robot.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();
  const double l_ref = 0.5 * diff_6d.transpose() * qi_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost->computeImpulseCost(robot, data, t, s), l_ref);
  cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_res);
  cost->computeImpulseCostHessian(robot, data, t, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::Jlog6(diff_SE3, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += J66_J_6d.transpose() * qi_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * qi_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  DerivativeChecker derivative_checker(robot);
  derivative_checker.setTestTolerance(1.0e-03);
  EXPECT_TRUE(derivative_checker.checkFirstOrderImpulseCostDerivatives(cost));
  // This is due to Gauss-Newton Hessian approximation.
  EXPECT_FALSE(derivative_checker.checkSecondOrderImpulseCostDerivatives(cost));
}


TEST_F(TaskSpace6DCostTest, fixedBase) {
  const int frame_id = 18;
  Robot robot(fixed_base_urdf);
  testStageCost(robot, frame_id);
  testTerminalCost(robot, frame_id);
  testImpulseCost(robot, frame_id);
}


TEST_F(TaskSpace6DCostTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf);
  for (const auto frame_id : frames) {
    testStageCost(robot, frame_id);
    testTerminalCost(robot, frame_id);
    testImpulseCost(robot, frame_id);
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}