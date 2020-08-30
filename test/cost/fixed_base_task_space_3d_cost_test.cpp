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

class FixedBaseTaskSpace3DCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    frame_id = 18;
    robot_ = Robot(urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    data_ = CostFunctionData(robot_);
    s = SplitSolution(robot_);
    kkt_res = KKTResidual(robot_);
    kkt_mat = KKTMatrix(robot_);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  int frame_id;
  std::string urdf_;
  Robot robot_;
  CostFunctionData data_;
  SplitSolution s;
  KKTResidual kkt_res;
  KKTMatrix kkt_mat;
};


TEST_F(FixedBaseTaskSpace3DCostTest, setWeights) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::Vector3d q_weight = Eigen::Vector3d::Random();
  const Eigen::Vector3d qf_weight = Eigen::Vector3d::Random();
  const Eigen::Vector3d q_ref = Eigen::Vector3d::Random();
  TaskSpace3DCost cost(robot_, frame_id);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_3d_weight(q_weight);
  cost.set_qf_3d_weight(qf_weight);
  cost.set_q_3d_ref(q_ref);
  s.q = Eigen::VectorXd::Random(dimq);
  s.v = Eigen::VectorXd::Random(dimv);
  s.a = Eigen::VectorXd::Random(dimv);
  s.f = Eigen::VectorXd::Random(robot_.max_dimf());
  s.u = Eigen::VectorXd::Random(dimv);
  robot_.updateKinematics(s.q, s.v, s.a);
  const Eigen::Vector3d q_task = robot_.framePosition(frame_id);
  const Eigen::Vector3d q_diff = q_task - q_ref;
  const double l_ref = dtau_ * 0.5 * q_diff.transpose() * q_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), l_ref);
  const double phi_ref = 0.5 * q_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  EXPECT_DOUBLE_EQ(cost.phi(robot_, data_, t_, s), phi_ref);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, robot_.dimv());
  robot_.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J_diff = robot_.frameRotation(frame_id) * J_6d.topRows(3);
  const Eigen::VectorXd lq_ref = dtau_ * J_diff.transpose() * q_weight.asDiagonal() * q_diff;
  cost.lq(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(lq_ref));
  const Eigen::VectorXd phiq_ref = J_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  kkt_res.lq().setZero();
  cost.phiq(robot_, data_, t_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(phiq_ref));
  const Eigen::MatrixXd lqq_ref = dtau_ * J_diff.transpose() * q_weight.asDiagonal() * J_diff;
  cost.lqq(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(lqq_ref));
  const Eigen::MatrixXd phiqq_ref = J_diff.transpose() * qf_weight.asDiagonal() * J_diff;
  kkt_mat.Qqq().setZero();
  cost.phiqq(robot_, data_, t_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(phiqq_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}