#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"


namespace idocp {

class StateConstraintRiccatiFactorizationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N;
};


void StateConstraintRiccatiFactorizationTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  StateConstraintRiccatiFactorization factorization(robot, N);
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), 0);
    EXPECT_EQ(factorization.T_aux(i).rows(), dimx);
    EXPECT_EQ(factorization.T_aux(i).cols(), 0);
  }
  EXPECT_EQ(factorization.ENEt().rows(), 0);
  EXPECT_EQ(factorization.ENEt().cols(), 0);
  EXPECT_EQ(factorization.EqNqq().rows(), 0);
  EXPECT_EQ(factorization.EqNqq().cols(), dimv);
  auto impulse_status = robot.createImpulseStatus();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    std::random_device rnd;
    if (rnd()%2==0) {
      impulse_status.activateImpulse(i);
    }
  }
  factorization.setImpulseStatus(impulse_status);
  const int dimf = impulse_status.dimp();
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), dimf);
    EXPECT_EQ(factorization.T_aux(i).rows(), dimx);
    EXPECT_EQ(factorization.T_aux(i).cols(), dimf);
  }
  EXPECT_EQ(factorization.ENEt().rows(), dimf);
  EXPECT_EQ(factorization.ENEt().cols(), dimf);
  EXPECT_EQ(factorization.EqNqq().rows(), dimf);
  EXPECT_EQ(factorization.EqNqq().cols(), dimv);
}


TEST_F(StateConstraintRiccatiFactorizationTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(StateConstraintRiccatiFactorizationTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}