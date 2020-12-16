#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class JointSpaceCostTest : public ::testing::Test {
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

  void test(Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau, t;
};


void JointSpaceCostTest::test(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
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
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimu);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimu);
  JointSpaceCost cost(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost.useKinematics());
  cost.set_q_weight(q_weight);
  cost.set_v_weight(v_weight);
  cost.set_a_weight(a_weight);
  cost.set_u_weight(u_weight);
  cost.set_qf_weight(qf_weight);
  cost.set_vf_weight(vf_weight);
  cost.set_q_ref(q_ref);
  cost.set_v_ref(v_ref);
  cost.set_a_ref(a_ref);
  cost.set_u_ref(u_ref);
  const SplitSolution s = SplitSolution::Random(robot);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double l_ref = 0.5 * dtau 
                           * ((q_weight.array()*q_diff.array()*q_diff.array()).sum()
                            + (v_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum()
                            + (a_weight.array()* (s.a-a_ref).array()*(s.a-a_ref).array()).sum()
                            + (u_weight.array()* (s.u-u_ref).array()*(s.u-u_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), l_ref);
  const double phi_ref = 0.5 * ((qf_weight.array()* q_diff.array()*q_diff.array()).sum()
                                + (vf_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.phi(robot, data, t, s), phi_ref);
  cost.lq(robot, data, t, dtau, s, kkt_res);
  cost.lv(robot, data, t, dtau, s, kkt_res);
  cost.la(robot, data, t, dtau, s, kkt_res);
  cost.lu(robot, data, t, dtau, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += dtau * Jq_diff.transpose() * q_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += dtau * q_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += dtau * v_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.la += dtau * a_weight.asDiagonal() * (s.a-a_ref);
  kkt_res_ref.lu() += dtau * u_weight.asDiagonal() * (s.u-u_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost.phiq(robot, data, t, s, kkt_res);
  cost.phiv(robot, data, t, s, kkt_res);
  if (robot.hasFloatingBase()) {
    kkt_res_ref.lq() += Jq_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += qf_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += vf_weight.asDiagonal() * (s.v-v_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost.lqq(robot, data, t, dtau, s, kkt_mat);
  cost.lvv(robot, data, t, dtau, s, kkt_mat);
  cost.laa(robot, data, t, dtau, s, kkt_mat);
  cost.luu(robot, data, t, dtau, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += dtau * Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += dtau * q_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += dtau * v_weight.asDiagonal();
  kkt_mat_ref.Qaa() += dtau * a_weight.asDiagonal();
  kkt_mat_ref.Quu() += dtau * u_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  cost.phiqq(robot, data, t, s, kkt_mat);
  cost.phivv(robot, data, t, s, kkt_mat);
  if (robot.hasFloatingBase()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * qf_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += qf_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += vf_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(JointSpaceCostTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  test(robot);
}


TEST_F(JointSpaceCostTest, floatingBase) {
  Robot robot(floating_base_urdf);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}