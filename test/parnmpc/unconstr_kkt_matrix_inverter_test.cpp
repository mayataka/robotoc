#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/parnmpc/unconstr_kkt_matrix_inverter.hpp"

#include "robot_factory.hpp"


namespace idocp {

class UnconstrKKTMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
};


TEST_F(UnconstrKKTMatrixInverterTest, test) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimH = 3*robot.dimv();
  const int dimKKT = 5*robot.dimv();
  const Eigen::MatrixXd H_seed_mat = Eigen::MatrixXd::Random(dimH, dimH);
  const Eigen::MatrixXd H_mat = H_seed_mat * H_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimH, dimH);
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  UnconstrKKTMatrixInverter inverter(robot);
  inverter.invert(dt, H_mat, KKT_mat_inv);
  Eigen::MatrixXd KKT_mat_ref = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  KKT_mat_ref.bottomRightCorner(dimH, dimH) = H_mat;
  KKT_mat_ref.block(   0, 3*dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  KKT_mat_ref.block(   0, 4*dimv, dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  KKT_mat_ref.block(dimv, 2*dimv, dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  KKT_mat_ref.block(dimv, 4*dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  KKT_mat_ref.block(dimx, 0, 3*dimv, 2*dimv) = KKT_mat_ref.block(0, dimx, 2*dimv, 3*dimv).transpose();
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat_ref.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat_ref).isIdentity());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}