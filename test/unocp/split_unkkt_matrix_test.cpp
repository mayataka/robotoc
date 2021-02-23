#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"


namespace idocp {

class SplitUnKKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    robot = Robot(urdf);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  std::string urdf;
  Robot robot;
  double dt;
};


TEST_F(SplitUnKKTMatrixTest, test) {
  SplitUnKKTMatrix matrix(robot);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimQ = 3*robot.dimv();
  const int dimKKT = 5*robot.dimv();
  const Eigen::MatrixXd Q_seed_mat = Eigen::MatrixXd::Random(dimQ, dimQ);
  const Eigen::MatrixXd Q_mat = Q_seed_mat * Q_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimQ, dimQ);
  matrix.Q = Q_mat;
  EXPECT_TRUE(Q_mat.topLeftCorner(dimv, dimv).isApprox(matrix.Qaa()));
  EXPECT_TRUE(Q_mat.block(0, dimv, dimv, dimv).isApprox(matrix.Qaq()));
  EXPECT_TRUE(Q_mat.block(0, 2*dimv, dimv, dimv).isApprox(matrix.Qav()));
  EXPECT_TRUE(Q_mat.bottomRightCorner(dimx, dimx).isApprox(matrix.Qxx()));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}