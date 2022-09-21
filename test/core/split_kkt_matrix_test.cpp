#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class SplitKKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ContactStatus& contact_status);
  static void test_isApprox(const Robot& robot, const ContactStatus& contact_status);

  double dt;
};


void SplitKKTMatrixTest::test(const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2 * robot.dimv();
  const int dimf = contact_status.dimf();
  EXPECT_EQ(kkt_mat.dimf(), dimf);

  EXPECT_EQ(kkt_mat.Fxx.rows(), dimx);
  EXPECT_EQ(kkt_mat.Fxx.cols(), dimx);
  EXPECT_EQ(kkt_mat.Fqq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fqq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fqv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fqv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fvq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fvq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fvv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fvv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fvu.rows(), dimv);
  EXPECT_EQ(kkt_mat.Fvu.cols(), dimu);

  EXPECT_EQ(kkt_mat.Qxx.rows(), dimx);
  EXPECT_EQ(kkt_mat.Qxx.cols(), dimx);
  EXPECT_EQ(kkt_mat.Qqq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qqv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qvq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qvq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qvv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qvv().cols(), dimv);

  EXPECT_EQ(kkt_mat.Qxu.rows(), dimx);
  EXPECT_EQ(kkt_mat.Qxu.cols(), dimu);
  EXPECT_EQ(kkt_mat.Qqu().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqu().cols(), dimu);
  EXPECT_EQ(kkt_mat.Qvu().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qvu().cols(), dimu);

  EXPECT_EQ(kkt_mat.Quu.rows(), dimu);
  EXPECT_EQ(kkt_mat.Quu.cols(), dimu);

  EXPECT_EQ(kkt_mat.Qaa.rows(), dimv);
  EXPECT_EQ(kkt_mat.Qaa.cols(), dimv);
  EXPECT_EQ(kkt_mat.Qff().rows(), dimf);
  EXPECT_EQ(kkt_mat.Qff().cols(), dimf);
  EXPECT_EQ(kkt_mat.Qqf().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqf().cols(), dimf);

  if (robot.hasFloatingBase()) {
    EXPECT_EQ(kkt_mat.Fqq_prev.rows(), dimv);
    EXPECT_EQ(kkt_mat.Fqq_prev.cols(), dimv);
  }
  else {
    EXPECT_EQ(kkt_mat.Fqq_prev.rows(), 0);
    EXPECT_EQ(kkt_mat.Fqq_prev.cols(), 0);
  }

  EXPECT_TRUE(kkt_mat.Fxx.isZero());
  EXPECT_TRUE(kkt_mat.Fvu.isZero());
  EXPECT_TRUE(kkt_mat.Qxx.isZero());
  EXPECT_TRUE(kkt_mat.Qxu.isZero());
  EXPECT_TRUE(kkt_mat.Quu.isZero());
  EXPECT_TRUE(kkt_mat.Qaa.isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  EXPECT_TRUE(kkt_mat.Qqf().isZero());
  EXPECT_TRUE(kkt_mat.Fqq_prev.isZero());

  kkt_mat.Fxx.setRandom();
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qxu.setRandom();
  EXPECT_TRUE(kkt_mat.Fxx.topLeftCorner(dimv, dimv).isApprox(kkt_mat.Fqq()));
  EXPECT_TRUE(kkt_mat.Fxx.topRightCorner(dimv, dimv).isApprox(kkt_mat.Fqv()));
  EXPECT_TRUE(kkt_mat.Fxx.bottomLeftCorner(dimv, dimv).isApprox(kkt_mat.Fvq()));
  EXPECT_TRUE(kkt_mat.Fxx.bottomRightCorner(dimv, dimv).isApprox(kkt_mat.Fvv()));
  EXPECT_TRUE(kkt_mat.Qxx.topLeftCorner(dimv, dimv).isApprox(kkt_mat.Qqq()));
  EXPECT_TRUE(kkt_mat.Qxx.topRightCorner(dimv, dimv).isApprox(kkt_mat.Qqv()));
  EXPECT_TRUE(kkt_mat.Qxx.bottomLeftCorner(dimv, dimv).isApprox(kkt_mat.Qvq()));
  EXPECT_TRUE(kkt_mat.Qxx.bottomRightCorner(dimv, dimv).isApprox(kkt_mat.Qvv()));
  EXPECT_TRUE(kkt_mat.Qxu.topRows(dimv).isApprox(kkt_mat.Qqu()));
  EXPECT_TRUE(kkt_mat.Qxu.bottomRows(dimv).isApprox(kkt_mat.Qvu()));

  EXPECT_NO_THROW(
    std::cout << kkt_mat << std::endl;
  );
}


void SplitKKTMatrixTest::test_isApprox(const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2 * robot.dimv();
  const int dimf = contact_status.dimf();

  kkt_mat.Fxx.setRandom();
  kkt_mat.Fvu.setRandom();
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qxu.setRandom();
  kkt_mat.Quu.setRandom();
  kkt_mat.Qff().setRandom();
  kkt_mat.Qqf().setRandom();

  auto kkt_mat_ref = kkt_mat;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fxx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fxx = kkt_mat.Fxx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fvu.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fvu = kkt_mat.Fvu;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxx = kkt_mat.Qxx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxu.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxu = kkt_mat.Qxu;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Quu.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Quu = kkt_mat.Quu;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  if (dimf > 0) {
    kkt_mat_ref.Qff().setRandom();
    EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qff() = kkt_mat.Qff();
    EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qqf().setRandom();
    EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qqf() = kkt_mat.Qqf();
    EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  }
  kkt_mat_ref.fx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.fx = kkt_mat.fx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qtt = Eigen::VectorXd::Random(1)[0];
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qtt = kkt_mat.Qtt;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.hx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.hx = kkt_mat.hx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.hu.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.hu = kkt_mat.hu;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


TEST_F(SplitKKTMatrixTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto contact_status = robot.createContactStatus();
  test(robot, contact_status);
  test_isApprox(robot, contact_status);
  contact_status.activateContact(0);
  test(robot, contact_status);
  test_isApprox(robot, contact_status);
}


TEST_F(SplitKKTMatrixTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  test(robot, contact_status);
  test_isApprox(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status);
  test_isApprox(robot, contact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}