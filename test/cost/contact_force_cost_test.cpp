#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class ContactForceCostTest : public ::testing::Test {
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

  void commonTest(Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau, t;
};


void ContactForceCostTest::commonTest(Robot& robot) const {
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Random());
    f_ref.push_back(Eigen::Vector3d::Random());
  }
  ContactForceCost cost(robot);
  CostFunctionData data(robot);
  EXPECT_FALSE(cost.useKinematics());
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  ContactStatus contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    contact_status.deactivateContact(i);
  }
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), 0);
  cost.lf(robot, data, t, dtau, s, kkt_res);
  cost.lff(robot, data, t, dtau, s, kkt_mat);
  EXPECT_TRUE(kkt_res.lf().isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  contact_status.setRandom();
  s.setRandom(robot, contact_status);
  double l_ref = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      l_ref += (f_weight[i].array() * (s.f[i].array()-f_ref[i].array()) 
                                    * (s.f[i].array()-f_ref[i].array())).sum();
    }
  }
  EXPECT_DOUBLE_EQ(cost.l(robot, data, t, dtau, s), 0.5*dtau*l_ref);
  kkt_res.setContactStatus(contact_status);
  kkt_mat.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  kkt_mat.Qff().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  cost.lf(robot, data, t, dtau, s, kkt_res);
  cost.lff(robot, data, t, dtau, s, kkt_mat);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).array()
          += dtau * f_weight[i].array() * (s.f[i].array()-f_ref[i].array());
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack) += dtau * f_weight[i];
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(ContactForceCostTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  commonTest(robot);
}


TEST_F(ContactForceCostTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  commonTest(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}