#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_gain.hpp"


namespace idocp {

class RiccatiGainTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


void RiccatiGainTest::test(const Robot& robot) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiGain gain(robot);
  EXPECT_EQ(gain.K.rows(), dimu);
  EXPECT_EQ(gain.K.cols(), 2*dimv);
  EXPECT_EQ(gain.k.size(), dimu);
  const Eigen::MatrixXd K_ref = Eigen::MatrixXd::Random(dimu, 2*dimv);
  const Eigen::VectorXd k_ref = Eigen::VectorXd::Random(dimu);
  gain.K = K_ref;
  gain.k = k_ref;
  EXPECT_TRUE(K_ref.isApprox(gain.K));
  EXPECT_TRUE(K_ref.leftCols(dimv).isApprox(gain.Kq()));
  EXPECT_TRUE(K_ref.rightCols(dimv).isApprox(gain.Kv()));
  EXPECT_TRUE(k_ref.isApprox(gain.k));
}


TEST_F(RiccatiGainTest, fixed_base) {
  test(fixed_base_robot);
}


TEST_F(RiccatiGainTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}