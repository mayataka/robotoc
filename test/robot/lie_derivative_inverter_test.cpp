#include <vector>

#include <gtest/gtest.h>
#include "Eigen/LU"

#include "robotoc/robot/lie_derivative_inverter.hpp"


namespace robotoc {

class LieDerivativeInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }
};


TEST_F(LieDerivativeInverterTest, test) {
  Eigen::MatrixXd lie_der(Eigen::MatrixXd::Zero(6, 6));
  lie_der.setRandom();
  lie_der.bottomLeftCorner(3, 3).setZero();
  const Eigen::MatrixXd inv_ref = lie_der.inverse();
  LieDerivativeInverter inverter;
  Eigen::MatrixXd inv(Eigen::MatrixXd::Zero(6, 6));
  inverter.computeLieDerivativeInverse(lie_der, inv);
  EXPECT_TRUE(inv.isApprox(inv_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}