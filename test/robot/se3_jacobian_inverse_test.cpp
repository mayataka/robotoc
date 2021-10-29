#include <vector>

#include <gtest/gtest.h>
#include "Eigen/LU"

#include "robotoc/robot/se3_jacobian_inverse.hpp"


namespace robotoc {

class SE3JacobianInverseTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }
};


TEST_F(SE3JacobianInverseTest, test) {
  Eigen::MatrixXd lie_der(Eigen::MatrixXd::Zero(6, 6));
  lie_der.setRandom();
  lie_der.bottomLeftCorner(3, 3).setZero();
  const Eigen::MatrixXd inv_ref = lie_der.inverse();
  SE3JacobianInverse inverter;
  Eigen::MatrixXd inv(Eigen::MatrixXd::Zero(6, 6));
  inverter.compute(lie_der, inv);
  EXPECT_TRUE(inv.isApprox(inv_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}