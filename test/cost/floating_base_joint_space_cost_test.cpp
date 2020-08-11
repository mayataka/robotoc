#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class FloatingBaseJointSpaceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    robot_ = Robot(urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    data_ = CostFunctionData(robot_);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::string urdf_;
  Robot robot_;
  CostFunctionData data_;
};


TEST_F(FloatingBaseJointSpaceCostTest, setWeights) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot_.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimv);
  JointSpaceCost cost(robot_);
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
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  robot_.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd f = Eigen::VectorXd::Random(robot_.dimf());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv);
  robot_.subtractConfiguration(q, q_ref, q_diff);
  const double l_ref = 0.5 * dtau_ 
                           * ((q_weight.array()*q_diff.array()*q_diff.array()).sum()
                            + (v_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum()
                            + (a_weight.array()* (a-a_ref).array()*(a-a_ref).array()).sum()
                            + (u_weight.array()* (u-u_ref).array()*(u-u_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, q, v, a, f, u), l_ref);
  const double phi_ref = 0.5 * ((qf_weight.array()* q_diff.array()*q_diff.array()).sum()
                                + (vf_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.phi(robot_, data_, t_, q, v), phi_ref);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv);
  cost.lq(robot_, data_, t_, dtau_, q, v, a, lq);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  robot_.dSubtractdConfigurationPlus(q, q_ref, Jq_diff);
  lq_ref = dtau_ * Jq_diff.transpose() * q_weight.asDiagonal() * q_diff;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  cost.lv(robot_, data_, t_, dtau_, q, v, a, lv);
  lv_ref = dtau_ * v_weight.asDiagonal() * (v-v_ref);
  EXPECT_TRUE(lv.isApprox(lv_ref));
  cost.la(robot_, data_, t_, dtau_, q, v, a, la);
  la_ref = dtau_ * a_weight.asDiagonal() * (a-a_ref);
  EXPECT_TRUE(la.isApprox(la_ref));
  cost.lu(robot_, data_, t_, dtau_, u, lu);
  lu_ref = dtau_ * u_weight.asDiagonal() * (u-u_ref);
  EXPECT_TRUE(lu.isApprox(lu_ref));
  cost.phiq(robot_, data_, t_, q, v, lq);
  lq_ref = Jq_diff.transpose() * qf_weight.asDiagonal() * q_diff;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  cost.phiv(robot_, data_, t_, q, v, lv);
  lv_ref = vf_weight.asDiagonal() * (v-v_ref);
  EXPECT_TRUE(lv.isApprox(lv_ref));
  Eigen::MatrixXd lqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd lvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd laa = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd lqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd lvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd laa_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd luu_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  lqq_ref = dtau_ * Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  lvv_ref = dtau_ * v_weight.asDiagonal();
  laa_ref = dtau_ * a_weight.asDiagonal();
  luu_ref = dtau_ * u_weight.asDiagonal();
  cost.lqq(robot_, data_, t_, dtau_, q, v, a, lqq);
  EXPECT_TRUE(lqq.isApprox(lqq_ref));
  cost.lvv(robot_, data_, t_, dtau_, q, v, a, lvv);
  EXPECT_TRUE(lvv.isApprox(lvv_ref));
  cost.laa(robot_, data_, t_, dtau_, q, v, a, laa);
  EXPECT_TRUE(laa.isApprox(laa_ref));
  cost.luu(robot_, data_, t_, dtau_, u, luu);
  EXPECT_TRUE(luu.isApprox(luu_ref));
  cost.augment_lqq(robot_, data_, t_, dtau_, q, v, a, lqq);
  EXPECT_TRUE(lqq.isApprox(2*lqq_ref));
  cost.augment_lvv(robot_, data_, t_, dtau_, q, v, a, lvv);
  EXPECT_TRUE(lvv.isApprox(2*lvv_ref));
  cost.augment_laa(robot_, data_, t_, dtau_, q, v, a, laa);
  EXPECT_TRUE(laa.isApprox(2*laa_ref));
  cost.augment_luu(robot_, data_, t_, dtau_, u, luu);
  EXPECT_TRUE(luu.isApprox(2*luu_ref));
  lqq.setZero();
  lvv.setZero();
  lqq_ref.setZero();
  lvv_ref.setZero();
  lqq_ref = Jq_diff.transpose() * qf_weight.asDiagonal() * Jq_diff;
  lvv_ref = vf_weight.asDiagonal();
  cost.phiqq(robot_, data_, t_, q, v, lqq);
  EXPECT_TRUE(lqq.isApprox(lqq_ref));
  cost.phivv(robot_, data_, t_, q, v, lvv);
  EXPECT_TRUE(lvv.isApprox(lvv_ref));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}