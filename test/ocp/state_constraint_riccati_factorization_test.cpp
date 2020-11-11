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
    max_num_impulse = 10;
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse;
};


void StateConstraintRiccatiFactorizationTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  StateConstraintRiccatiFactorization factorization(robot, N, max_num_impulse);
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), 0);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_EQ(factorization.T_impulse(i).rows(), dimx);
    EXPECT_EQ(factorization.T_impulse(i).cols(), 0);
    EXPECT_EQ(factorization.T_lift(i).rows(), dimx);
    EXPECT_EQ(factorization.T_lift(i).cols(), 0);
  }
  EXPECT_EQ(factorization.Eq().rows(), 0);
  EXPECT_EQ(factorization.Eq().cols(), dimv);
  EXPECT_EQ(factorization.ENEt().rows(), 0);
  EXPECT_EQ(factorization.ENEt().cols(), 0);
  EXPECT_EQ(factorization.EN().rows(), 0);
  EXPECT_EQ(factorization.EN().cols(), dimx);
  EXPECT_EQ(factorization.ENq().rows(), 0);
  EXPECT_EQ(factorization.ENq().cols(), dimv);
  EXPECT_EQ(factorization.e().size(), 0);
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  factorization.setImpulseStatus(impulse_status);
  const int dimf = impulse_status.dimp();
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), dimf);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_EQ(factorization.T_impulse(i).rows(), dimx);
    EXPECT_EQ(factorization.T_impulse(i).cols(), dimf);
    EXPECT_EQ(factorization.T_lift(i).rows(), dimx);
    EXPECT_EQ(factorization.T_lift(i).cols(), dimf);
  }
  EXPECT_EQ(factorization.Eq().rows(), dimf);
  EXPECT_EQ(factorization.Eq().cols(), dimv);
  EXPECT_EQ(factorization.ENEt().rows(), dimf);
  EXPECT_EQ(factorization.ENEt().cols(), dimf);
  EXPECT_EQ(factorization.EN().rows(), dimf);
  EXPECT_EQ(factorization.EN().cols(), dimx);
  EXPECT_EQ(factorization.ENq().rows(), dimf);
  EXPECT_EQ(factorization.ENq().cols(), dimv);
  EXPECT_EQ(factorization.e().size(), dimf);
  factorization.EN().setRandom();
  EXPECT_TRUE(factorization.EN().leftCols(dimv).isApprox(factorization.ENq()));
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