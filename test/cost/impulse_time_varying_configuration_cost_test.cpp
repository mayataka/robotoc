#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

class ImpulseTimeVaryingConfigurationCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test(Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double t;
};


void ImpulseTimeVaryingConfigurationCostTest::test(Robot& robot) const {
  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  ImpulseSplitKKTMatrix kkt_mat(robot);
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_mat.Qqq().setRandom();
  kkt_mat.Qvv().setRandom();
  kkt_mat.Qdvdv().setRandom();
  kkt_res.lx().setRandom();
  kkt_res.ldv.setRandom();
  ImpulseSplitKKTMatrix kkt_mat_ref = kkt_mat;
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd dv_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd q0 = Eigen::VectorXd::Zero(dimq);
  robot.generateFeasibleConfiguration(q0);
  const Eigen::VectorXd v0 = Eigen::VectorXd::Random(dimv); 
  const double t0 = std::abs(Eigen::VectorXd::Random(1)[0]);
  ImpulseTimeVaryingConfigurationCost cost(robot);
  CostFunctionData data(robot);
  cost.set_q_weight(q_weight);
  cost.set_v_weight(v_weight);
  cost.set_dv_weight(dv_weight);
  cost.set_ref(t0, q0, v0);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  robot.integrateConfiguration(q0, v0, t-t0, q_ref);
  const Eigen::VectorXd v_ref = v0;
  const Eigen::VectorXd dv_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd q_diff = Eigen::VectorXd::Zero(dimv); 
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref, q_diff);
  }
  else {
    q_diff = s.q - q_ref;
  }
  const double l_ref = 0.5 * ((q_weight.array()*q_diff.array()*q_diff.array()).sum()
                            + (v_weight.array()* (s.v-v_ref).array()*(s.v-v_ref).array()).sum()
                            + (dv_weight.array()* (s.dv-dv_ref).array()*(s.dv-dv_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, s), l_ref);
  cost.lq(robot, data, t, s, kkt_res);
  cost.lv(robot, data, t, s, kkt_res);
  cost.ldv(robot, data, t, s, kkt_res);
  Eigen::MatrixXd Jq_diff = Eigen::MatrixXd::Zero(dimv, dimv);
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref, Jq_diff);
    kkt_res_ref.lq() += Jq_diff.transpose() * q_weight.asDiagonal() * q_diff;
  }
  else {
    kkt_res_ref.lq() += q_weight.asDiagonal() * (s.q-q_ref);
  }
  kkt_res_ref.lv() += v_weight.asDiagonal() * (s.v-v_ref);
  kkt_res_ref.ldv += dv_weight.asDiagonal() * (s.dv-dv_ref);
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  cost.lqq(robot, data, t, s, kkt_mat);
  cost.lvv(robot, data, t, s, kkt_mat);
  cost.ldvdv(robot, data, t, s, kkt_mat);
  if (robot.has_floating_base()) {
    kkt_mat_ref.Qqq() += Jq_diff.transpose() * q_weight.asDiagonal() * Jq_diff;
  }
  else {
    kkt_mat_ref.Qqq() += q_weight.asDiagonal();
  }
  kkt_mat_ref.Qvv() += v_weight.asDiagonal();
  kkt_mat_ref.Qdvdv() += dv_weight.asDiagonal();
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(ImpulseTimeVaryingConfigurationCostTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  test(robot);
}


TEST_F(ImpulseTimeVaryingConfigurationCostTest, floatingBase) {
  Robot robot(floating_base_urdf);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}