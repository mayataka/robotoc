#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/task_space_6d_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


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

  void setRefByTransAndRot(Robot& robot, const int frame_id) const;
  void setRefBySE3(Robot& robot, const int frame_id) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau, t;
};


void TaskSpace6DCostTest::setRefByTransAndRot(Robot& robot, const int frame_id) const {
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
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  const Eigen::Vector3d position_ref = ref_placement.translation();
  const Eigen::Matrix3d rotation_ref = ref_placement.rotation();
  TaskSpace6DCost cost(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost.set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost.set_q_6d_ref(position_ref, rotation_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();

  const double l_ref = dtau * 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), l_ref);
  const double phi_ref = 0.5 * diff_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.phi(robot, data, t, s), phi_ref);

  cost.lq(robot, data, t, dtau, s, kkt_res);
  cost.lqq(robot, data, t, dtau, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::Jlog6(diff_SE3, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost.phiq(robot, data, t, s, kkt_res);
  cost.phiqq(robot, data, t, s, kkt_mat);
  kkt_res_ref.lq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void TaskSpace6DCostTest::setRefBySE3(Robot& robot, const int frame_id) const {
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
  const pinocchio::SE3 ref_placement = pinocchio::SE3::Random();
  TaskSpace6DCost cost(robot, frame_id);
  CostFunctionData data(robot);
  EXPECT_TRUE(cost.useKinematics());
  cost.set_q_6d_weight(q_weight.tail(3), q_weight.head(3));
  cost.set_qf_6d_weight(qf_weight.tail(3), qf_weight.head(3));
  cost.set_q_6d_ref(ref_placement);
  const SplitSolution s = SplitSolution::Random(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const pinocchio::SE3 placement = robot.framePlacement(frame_id);
  const pinocchio::SE3 diff_SE3 = ref_placement.inverse() * placement;
  const Eigen::VectorXd diff_6d = pinocchio::log6(diff_SE3).toVector();

  const double l_ref = dtau * 0.5 * diff_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), l_ref);
  const double phi_ref = 0.5 * diff_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  EXPECT_DOUBLE_EQ(cost.phi(robot, data, t, s), phi_ref);

  cost.lq(robot, data, t, dtau, s, kkt_res);
  cost.lqq(robot, data, t, dtau, s, kkt_mat);
  Eigen::MatrixXd J_66 = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd J_6d = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::Jlog6(diff_SE3, J_66);
  robot.getFrameJacobian(frame_id, J_6d);
  const Eigen::MatrixXd J66_J_6d = J_66 * J_6d;
  kkt_res_ref.lq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += dtau * J66_J_6d.transpose() * q_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost.phiq(robot, data, t, s, kkt_res);
  cost.phiqq(robot, data, t, s, kkt_mat);
  kkt_res_ref.lq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * diff_6d;
  kkt_mat_ref.Qqq() += J66_J_6d.transpose() * qf_weight.asDiagonal() * J66_J_6d;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(TaskSpace6DCostTest, fixedBase) {
  const int frame_id = 18;
  Robot robot(fixed_base_urdf);
  setRefByTransAndRot(robot, frame_id);
  setRefBySE3(robot, frame_id);
}


TEST_F(TaskSpace6DCostTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf);
  for (const auto frame_id : frames) {
    setRefByTransAndRot(robot, frame_id);
    setRefBySE3(robot, frame_id);
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}