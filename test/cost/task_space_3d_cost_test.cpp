#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class TaskSpace3DCostTest : public ::testing::Test {
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

  void test(Robot& robot, const int frame_id) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau, t;
};


void TaskSpace3DCostTest::test(Robot& robot, const int frame_id) const {
  const int dimv = robot.dimv();
  KKTMatrix kkt_mat(robot);
  KKTResidual kkt_res(robot);
  kkt_mat.Qqq().setRandom();
  kkt_mat.Qvv().setRandom();
  kkt_mat.Qaa().setRandom();
  kkt_mat.Quu().setRandom();
  kkt_res.lq().setRandom();
  kkt_res.lv().setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu().setRandom();
  KKTMatrix kkt_mat_ref = kkt_mat;
  KKTResidual kkt_res_ref = kkt_res;
  const Eigen::Vector3d q_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d q_ref = Eigen::Vector3d::Random();
  TaskSpace3DCost cost(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_3d_weight(q_weight);
  cost.set_qf_3d_weight(qf_weight);
  cost.set_q_3d_ref(q_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_task = robot.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - q_ref;
  const double l_ref = dtau * 0.5 * q_diff.transpose() * q_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), l_ref);
  const double phi_ref = 0.5 * q_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost.phi(robot, data, t, s), phi_ref);

  cost.lq(robot, data, t, dtau, s, kkt_res);
  cost.lqq(robot, data, t, dtau, s, kkt_mat);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot.frameRotation(frame_id) * J_6d.topRows(3);
  kkt_res_ref.lq() += dtau * J_diff.transpose() * q_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += dtau * J_diff.transpose() * q_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost.phiq(robot, data, t, s, kkt_res);
  cost.phiqq(robot, data, t, s, kkt_mat);
  kkt_res_ref.lq() += J_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  kkt_mat_ref.Qqq() += J_diff.transpose() * qf_weight.asDiagonal() * J_diff;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(TaskSpace3DCostTest, fixedBase) {
  const int frame_id = 18;
  Robot robot(fixed_base_urdf);
  test(robot, frame_id);
}


TEST_F(TaskSpace3DCostTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf);
  for (const auto frame_id : frames) {
    test(robot, frame_id);
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}