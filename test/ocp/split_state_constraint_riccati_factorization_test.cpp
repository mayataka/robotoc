#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_state_constraint_riccati_factorization.hpp"


namespace idocp {

class SplitStateConstraintRiccatiFactorizationTest : public ::testing::Test {
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

  void testDim(const Robot& robot) const;
  void testAssignment(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse;
};


void SplitStateConstraintRiccatiFactorizationTest::testDim(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  SplitStateConstraintRiccatiFactorization factorization(robot, N, max_num_impulse);
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), 0);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_EQ(factorization.T_impulse(i).rows(), dimx);
    EXPECT_EQ(factorization.T_impulse(i).cols(), 0);
    EXPECT_EQ(factorization.T_aux(i).rows(), dimx);
    EXPECT_EQ(factorization.T_aux(i).cols(), 0);
    EXPECT_EQ(factorization.T_lift(i).rows(), dimx);
    EXPECT_EQ(factorization.T_lift(i).cols(), 0);
  }
  EXPECT_EQ(factorization.Eq().rows(), 0);
  EXPECT_EQ(factorization.Eq().cols(), dimv);
  EXPECT_EQ(factorization.EN().rows(), 0);
  EXPECT_EQ(factorization.EN().cols(), dimx);
  EXPECT_EQ(factorization.ENq().rows(), 0);
  EXPECT_EQ(factorization.ENq().cols(), dimv);
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
    EXPECT_EQ(factorization.T_aux(i).rows(), dimx);
    EXPECT_EQ(factorization.T_aux(i).cols(), dimf);
    EXPECT_EQ(factorization.T_lift(i).rows(), dimx);
    EXPECT_EQ(factorization.T_lift(i).cols(), dimf);
  }
  EXPECT_EQ(factorization.Eq().rows(), dimf);
  EXPECT_EQ(factorization.Eq().cols(), dimv);
  EXPECT_EQ(factorization.EN().rows(), dimf);
  EXPECT_EQ(factorization.EN().cols(), dimx);
  EXPECT_EQ(factorization.ENq().rows(), dimf);
  EXPECT_EQ(factorization.ENq().cols(), dimv);
  factorization.EN().setRandom();
  EXPECT_TRUE(factorization.EN().leftCols(dimv).isApprox(factorization.ENq()));
  factorization.resetImpulseStatus();
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(factorization.T(i).rows(), dimx);
    EXPECT_EQ(factorization.T(i).cols(), 0);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_EQ(factorization.T_impulse(i).rows(), dimx);
    EXPECT_EQ(factorization.T_impulse(i).cols(), 0);
    EXPECT_EQ(factorization.T_aux(i).rows(), dimx);
    EXPECT_EQ(factorization.T_aux(i).cols(), 0);
    EXPECT_EQ(factorization.T_lift(i).rows(), dimx);
    EXPECT_EQ(factorization.T_lift(i).cols(), 0);
  }
  EXPECT_EQ(factorization.Eq().rows(), 0);
  EXPECT_EQ(factorization.Eq().cols(), dimv);
  EXPECT_EQ(factorization.EN().rows(), 0);
  EXPECT_EQ(factorization.EN().cols(), dimx);
  EXPECT_EQ(factorization.ENq().rows(), 0);
  EXPECT_EQ(factorization.ENq().cols(), dimv);
}


void SplitStateConstraintRiccatiFactorizationTest::testAssignment(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  SplitStateConstraintRiccatiFactorization factorization(robot, N, max_num_impulse);
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  factorization.setImpulseStatus(impulse_status);
  const int dimf = impulse_status.dimp();
  std::vector<Eigen::MatrixXd> T;
  for (int i=0; i<N; ++i) {
    T.push_back(Eigen::MatrixXd::Random(dimx, dimf));
    factorization.T(i) = T[i];
  }
  std::vector<Eigen::MatrixXd> T_impulse, T_aux, T_lift;
  for (int i=0; i<max_num_impulse; ++i) {
    T_impulse.push_back(Eigen::MatrixXd::Random(dimx, dimf));
    T_aux.push_back(Eigen::MatrixXd::Random(dimx, dimf));
    T_lift.push_back(Eigen::MatrixXd::Random(dimx, dimf));
    factorization.T_impulse(i) = T_impulse[i];
    factorization.T_aux(i) = T_aux[i];
    factorization.T_lift(i) = T_lift[i];
  }
  const Eigen::MatrixXd Eq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd EN = Eigen::MatrixXd::Random(dimf, dimx);
  factorization.Eq() = Eq;
  factorization.EN() = EN;
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(factorization.T(i).isApprox(T[i]));
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.T_impulse(i).isApprox(T_impulse[i]));
    EXPECT_TRUE(factorization.T_aux(i).isApprox(T_aux[i]));
    EXPECT_TRUE(factorization.T_lift(i).isApprox(T_lift[i]));
  }
  EXPECT_TRUE(factorization.Eq().isApprox(Eq));
  EXPECT_TRUE(factorization.EN().isApprox(EN));
  EXPECT_TRUE(factorization.ENq().isApprox(EN.leftCols(dimv)));
}


TEST_F(SplitStateConstraintRiccatiFactorizationTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testDim(robot);
  testAssignment(robot);
}


TEST_F(SplitStateConstraintRiccatiFactorizationTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testDim(robot);
  testAssignment(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}