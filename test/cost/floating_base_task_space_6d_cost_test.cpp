#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/task_space_6d_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class FloatingBaseTaskSpace6DCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    frame_id = 14;
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


TEST_F(FloatingBaseTaskSpace6DCostTest, setRefByVectorAndMatrix) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(6);
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  TaskSpace6DCost cost(robot_, frame_id);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost.set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost.set_q_6d_ref(position_ref, rotation_ref);
  s.q = Eigen::VectorXd::Random(dimq);
  s.v = Eigen::VectorXd::Random(dimv);
  s.a = Eigen::VectorXd::Random(dimv);
  s.f = Eigen::VectorXd::Random(robot_.max_dimf());
  s.u = Eigen::VectorXd::Random(dimv);
  robot_.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot_.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();
  const double l_ref = dtau_ * 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), l_ref);
  const double phi_ref = 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.phi(robot_, data_, t_, s), phi_ref);

  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, robot_.dimv());
  pinocchio::Jlog6(diff_SE3, J_66);
  robot_.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  Eigen::VectorXd lq = dtau_ * J66_J_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  cost.lq(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(lq));
  Eigen::MatrixXd lqq = dtau_ * J66_J_6d.transpose() * q_weight.asDiagonal() * J66_J_6d;
  cost.lqq(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(lqq));

  lq += J66_J_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  cost.phiq(robot_, data_, t_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(lq));
  lqq += J66_J_6d.transpose() * qf_weight.asDiagonal() * J66_J_6d;
  cost.phiqq(robot_, data_, t_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(lqq));
}


TEST_F(FloatingBaseTaskSpace6DCostTest, setRefBySE3) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(6);
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  TaskSpace6DCost cost(robot_, frame_id);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost.set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost.set_q_6d_ref(ref_placement);
  s.q = Eigen::VectorXd::Random(dimq);
  s.v = Eigen::VectorXd::Random(dimv);
  s.a = Eigen::VectorXd::Random(dimv);
  s.f = Eigen::VectorXd::Random(robot_.max_dimf());
  s.u = Eigen::VectorXd::Random(dimv);
  robot_.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot_.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();
  const double l_ref = dtau_ * 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, s), l_ref);
  const double phi_ref = 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.phi(robot_, data_, t_, s), phi_ref);

  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, robot_.dimv());
  pinocchio::Jlog6(diff_SE3, J_66);
  robot_.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  Eigen::VectorXd lq = dtau_ * J66_J_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  cost.lq(robot_, data_, t_, dtau_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(lq));
  Eigen::MatrixXd lqq = dtau_ * J66_J_6d.transpose() * q_weight.asDiagonal() * J66_J_6d;
  cost.lqq(robot_, data_, t_, dtau_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(lqq));

  lq += J66_J_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  cost.phiq(robot_, data_, t_, s, kkt_res);
  EXPECT_TRUE(kkt_res.lq().isApprox(lq));
  lqq += J66_J_6d.transpose() * qf_weight.asDiagonal() * J66_J_6d;
  cost.phiqq(robot_, data_, t_, s, kkt_mat);
  EXPECT_TRUE(kkt_mat.Qqq().isApprox(lqq));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}