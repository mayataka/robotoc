#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"


namespace idocp {

class SplitKKTMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;
  void testWithImpulse(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double dtau;
};


void SplitKKTMatrixInverterTest::test(const Robot& robot) const {
  SplitKKTMatrix matrix(robot);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  const int dimQ = 2*robot.dimv() + robot.dimu();
  const int dimKKT = 4*robot.dimv() + robot.dimu();
  const Eigen::MatrixXd KKT_seed_mat = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd KKT_mat = KKT_seed_mat * KKT_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  KKT_mat.topLeftCorner(dimx, dimx).setZero();
  KKT_mat.topRows(dimv).setZero();
  auto q_prev = robot.generateFeasibleConfiguration();
  auto q_next = robot.generateFeasibleConfiguration();
  robot.dSubtractdConfigurationMinus(q_prev, q_next, KKT_mat.block(0, dimx+dimu, dimv, dimv));
  KKT_mat.block(             0, dimx+dimu+dimv, dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    KKT_mat.block(             0, dimx+dimu+dimv, 6, 6).setRandom();
  }
  KKT_mat.bottomLeftCorner(dimQ, dimx) = KKT_mat.topRightCorner(dimx, dimQ).transpose();
  matrix.Fqu() = KKT_mat.block(             0,           dimx, dimv, dimu);
  matrix.Fqq() = KKT_mat.block(             0,      dimx+dimu, dimv, dimv);
  matrix.Fqv() = KKT_mat.block(             0, dimx+dimu+dimv, dimv, dimv);
  matrix.Fvu() = KKT_mat.block(          dimv,           dimx, dimv, dimu);
  matrix.Fvq() = KKT_mat.block(          dimv,      dimx+dimu, dimv, dimv);
  matrix.Fvv() = KKT_mat.block(          dimv, dimx+dimu+dimv, dimv, dimv);
  matrix.Quu() = KKT_mat.block(          dimx,           dimx, dimu, dimu);
  matrix.Quq() = KKT_mat.block(          dimx,      dimx+dimu, dimu, dimv);
  matrix.Quv() = KKT_mat.block(          dimx, dimx+dimu+dimv, dimu, dimv);
  matrix.Qqu() = KKT_mat.block(     dimx+dimu,           dimx, dimv, dimu);
  matrix.Qqq() = KKT_mat.block(     dimx+dimu,      dimx+dimu, dimv, dimv);
  matrix.Qqv() = KKT_mat.block(     dimx+dimu, dimx+dimu+dimv, dimv, dimv);
  matrix.Qvu() = KKT_mat.block(dimx+dimu+dimv,           dimx, dimv, dimu);
  matrix.Qvq() = KKT_mat.block(dimx+dimu+dimv,      dimx+dimu, dimv, dimv);
  matrix.Qvv() = KKT_mat.block(dimx+dimu+dimv, dimx+dimu+dimv, dimv, dimv);
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(dtau, matrix.Jac(), matrix.Qss(), KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity());
}


void SplitKKTMatrixInverterTest::testWithImpulse(const Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  if (robot.hasFloatingBase()) {
    std::random_device rnd;
    for (int i=0; i<robot.maxPointContacts(); ++i) {
      if (rnd() % 2 == 0) {
        impulse_status.activateImpulse(i);
      }
    }
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  else {
    impulse_status.activateImpulse(0);
  }
  SplitKKTMatrix matrix(robot);
  matrix.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  const int dimQ = 2*robot.dimv() + robot.dimu();
  const int dimi = impulse_status.dimf();
  const int dimKKT = 4*robot.dimv() + robot.dimu() + dimi;
  ASSERT_TRUE(dimi > 0);
  const Eigen::MatrixXd KKT_seed_mat = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd KKT_mat = KKT_seed_mat * KKT_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  KKT_mat.topLeftCorner(dimx+dimi, dimx+dimi).setZero();
  KKT_mat.topRows(dimv).setZero();
  KKT_mat.middleRows(dimx, dimi).setZero();
  auto q_prev = robot.generateFeasibleConfiguration();
  auto q_next = robot.generateFeasibleConfiguration();
  robot.dSubtractdConfigurationMinus(q_prev, q_next, KKT_mat.block(0, dimx+dimi+dimu, dimv, dimv));
  KKT_mat.block(             0, dimx+dimi+dimu+dimv, dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    KKT_mat.block(             0, dimx+dimi+dimu+dimv, 6, 6).setRandom();
  }
  KKT_mat.block(          dimx,      dimx+dimi+dimu, dimi, dimv).setRandom();
  KKT_mat.bottomLeftCorner(dimQ, dimx+dimi) = KKT_mat.topRightCorner(dimx+dimi, dimQ).transpose();
  matrix.Fqu() = KKT_mat.block(                  0,           dimx+dimi, dimv, dimu);
  matrix.Fqq() = KKT_mat.block(                  0,      dimx+dimi+dimu, dimv, dimv);
  matrix.Fqv() = KKT_mat.block(                  0, dimx+dimi+dimu+dimv, dimv, dimv);
  matrix.Fvu() = KKT_mat.block(               dimv,           dimx+dimi, dimv, dimu);
  matrix.Fvq() = KKT_mat.block(               dimv,      dimx+dimi+dimu, dimv, dimv);
  matrix.Fvv() = KKT_mat.block(               dimv, dimx+dimi+dimu+dimv, dimv, dimv);
  matrix.Pq()  = KKT_mat.block(               dimx,      dimx+dimi+dimu, dimi, dimv);
  matrix.Quu() = KKT_mat.block(          dimx+dimi,           dimx+dimi, dimu, dimu);
  matrix.Quq() = KKT_mat.block(          dimx+dimi,      dimx+dimi+dimu, dimu, dimv);
  matrix.Quv() = KKT_mat.block(          dimx+dimi, dimx+dimi+dimu+dimv, dimu, dimv);
  matrix.Qqu() = KKT_mat.block(     dimx+dimi+dimu,           dimx+dimi, dimv, dimu);
  matrix.Qqq() = KKT_mat.block(     dimx+dimi+dimu,      dimx+dimi+dimu, dimv, dimv);
  matrix.Qqv() = KKT_mat.block(     dimx+dimi+dimu, dimx+dimi+dimu+dimv, dimv, dimv);
  matrix.Qvu() = KKT_mat.block(dimx+dimi+dimu+dimv,           dimx+dimi, dimv, dimu);
  matrix.Qvq() = KKT_mat.block(dimx+dimi+dimu+dimv,      dimx+dimi+dimu, dimv, dimv);
  matrix.Qvv() = KKT_mat.block(dimx+dimi+dimu+dimv, dimx+dimi+dimu+dimv, dimv, dimv);
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(dtau, matrix.Jac(), matrix.Pq(), matrix.Qss(), KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref, 1.0e-08));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity(1.0e-06));
}


TEST_F(SplitKKTMatrixInverterTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  test(robot);
  testWithImpulse(robot);
}


TEST_F(SplitKKTMatrixInverterTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
  testWithImpulse(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}