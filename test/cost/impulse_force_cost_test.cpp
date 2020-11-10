#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

class ImpulseForceCostTest : public ::testing::Test {
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

  void commonTest(Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double t;
};


void ImpulseForceCostTest::commonTest(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
  }
  ImpulseForceCost cost(robot);
  CostFunctionData data(robot);
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    impulse_status.deactivateImpulse(i);
  }
  ImpulseKKTMatrix kkt_mat(robot);
  ImpulseKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_mat.setImpulseStatus(impulse_status);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, s), 0);
  cost.lf(robot, data, t, s, kkt_res);
  cost.lff(robot, data, t, s, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  std::random_device rnd;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (rnd()%2 == 0) 
      impulse_status.activateImpulse(i);
  }
  s.setRandom(robot, impulse_status);
  double l_ref = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      l_ref += (f_weight[i].array() * (s.f[i].array()-f_ref[i].array()) 
                                    * (s.f[i].array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, s), 0.5*l_ref);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_mat.setImpulseStatus(impulse_status);
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  ImpulseKKTResidual kkt_res_ref = kkt_res;
  ImpulseKKTMatrix kkt_mat_ref = kkt_mat;
  cost.lf(robot, data, t, s, kkt_res);
  cost.lff(robot, data, t, s, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += f_weight[i].array() * (s.f[i].array()-f_ref[i].array());
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += f_weight[i];
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(ImpulseForceCostTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  commonTest(robot);
}


TEST_F(ImpulseForceCostTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  commonTest(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}