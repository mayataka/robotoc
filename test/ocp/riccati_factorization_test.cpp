#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

class RiccatiFactorizationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
    fixed_base_robot_contact = Robot(fixed_base_urdf, {18});
    floating_base_robot_contact = Robot(floating_base_urdf, {14, 24, 34, 44});
  }

  virtual void TearDown() {
  }

  static void testWithoutContact(const Robot& robot);
  static void testWithContact(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
  Robot fixed_base_robot_contact, floating_base_robot_contact;
};


void RiccatiFactorizationTest::testWithoutContact(const Robot& robot) {
  ASSERT_TRUE(robot.max_point_contacts() == 0);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  RiccatiFactorization riccati(robot);
  EXPECT_EQ(riccati.Pqq.rows(), dimv);
  EXPECT_EQ(riccati.Pqq.cols(), dimv);
  EXPECT_EQ(riccati.Pqv.rows(), dimv);
  EXPECT_EQ(riccati.Pqv.cols(), dimv);
  EXPECT_EQ(riccati.Pvq.rows(), dimv);
  EXPECT_EQ(riccati.Pvq.cols(), dimv);
  EXPECT_EQ(riccati.Pvv.rows(), dimv);
  EXPECT_EQ(riccati.Pvv.cols(), dimv);
  EXPECT_EQ(riccati.sq.size(), dimv);
  EXPECT_EQ(riccati.sv.size(), dimv);
  EXPECT_EQ(riccati.K.rows(), dimu);
  EXPECT_EQ(riccati.K.cols(), dimx);
  EXPECT_EQ(riccati.k.size(), dimu);
  EXPECT_EQ(riccati.ApBK.rows(), 0);
  EXPECT_EQ(riccati.ApBK.cols(), 0);
  EXPECT_EQ(riccati.apBk.size(), 0);
  EXPECT_EQ(riccati.BGinvBt.rows(), 0);
  EXPECT_EQ(riccati.BGinvBt.cols(), 0);
  EXPECT_EQ(riccati.Pi.rows(), 0);
  EXPECT_EQ(riccati.Pi.cols(), 0);
  EXPECT_EQ(riccati.pi.size(), 0);
  EXPECT_EQ(riccati.N.rows(), 0);
  EXPECT_EQ(riccati.N.cols(), 0);
  EXPECT_EQ(riccati.Kq().rows(), dimu);
  EXPECT_EQ(riccati.Kq().cols(), dimv);
  EXPECT_EQ(riccati.Kv().rows(), dimu);
  EXPECT_EQ(riccati.Kv().cols(), dimv);
  const Eigen::MatrixXd K_ref = Eigen::MatrixXd::Random(dimu, dimx);
  const Eigen::VectorXd k_ref = Eigen::VectorXd::Random(dimu);
  riccati.K = K_ref;
  riccati.k = k_ref;
  EXPECT_TRUE(K_ref.isApprox(riccati.K));
  EXPECT_TRUE(k_ref.isApprox(riccati.k));
  EXPECT_TRUE(K_ref.leftCols(dimv).isApprox(riccati.Kq()));
  EXPECT_TRUE(K_ref.rightCols(dimv).isApprox(riccati.Kv()));
}


void RiccatiFactorizationTest::testWithContact(const Robot& robot) {
  ASSERT_TRUE(robot.max_point_contacts() > 0);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  RiccatiFactorization riccati(robot);
  EXPECT_EQ(riccati.Pqq.rows(), dimv);
  EXPECT_EQ(riccati.Pqq.cols(), dimv);
  EXPECT_EQ(riccati.Pqv.rows(), dimv);
  EXPECT_EQ(riccati.Pqv.cols(), dimv);
  EXPECT_EQ(riccati.Pvq.rows(), dimv);
  EXPECT_EQ(riccati.Pvq.cols(), dimv);
  EXPECT_EQ(riccati.Pvv.rows(), dimv);
  EXPECT_EQ(riccati.Pvv.cols(), dimv);
  EXPECT_EQ(riccati.sq.size(), dimv);
  EXPECT_EQ(riccati.sv.size(), dimv);
  EXPECT_EQ(riccati.K.rows(), dimu);
  EXPECT_EQ(riccati.K.cols(), dimx);
  EXPECT_EQ(riccati.k.size(), dimu);
  EXPECT_EQ(riccati.ApBK.rows(), dimx);
  EXPECT_EQ(riccati.ApBK.cols(), dimx);
  EXPECT_EQ(riccati.apBk.size(), dimx);
  EXPECT_EQ(riccati.BGinvBt.rows(), dimv);
  EXPECT_EQ(riccati.BGinvBt.cols(), dimv);
  EXPECT_EQ(riccati.Pi.rows(), dimx);
  EXPECT_EQ(riccati.Pi.cols(), dimx);
  EXPECT_EQ(riccati.pi.size(), dimx);
  EXPECT_EQ(riccati.N.rows(), dimx);
  EXPECT_EQ(riccati.N.cols(), dimx);
  EXPECT_EQ(riccati.Kq().rows(), dimu);
  EXPECT_EQ(riccati.Kq().cols(), dimv);
  EXPECT_EQ(riccati.Kv().rows(), dimu);
  EXPECT_EQ(riccati.Kv().cols(), dimv);
  const Eigen::MatrixXd K_ref = Eigen::MatrixXd::Random(dimu, dimx);
  const Eigen::VectorXd k_ref = Eigen::VectorXd::Random(dimu);
  riccati.K = K_ref;
  riccati.k = k_ref;
  EXPECT_TRUE(K_ref.isApprox(riccati.K));
  EXPECT_TRUE(k_ref.isApprox(riccati.k));
  EXPECT_TRUE(K_ref.leftCols(dimv).isApprox(riccati.Kq()));
  EXPECT_TRUE(K_ref.rightCols(dimv).isApprox(riccati.Kv()));
}


TEST_F(RiccatiFactorizationTest, fixed_base) {
  testWithoutContact(fixed_base_robot);
  testWithContact(fixed_base_robot_contact);
}


TEST_F(RiccatiFactorizationTest, floating_base) {
  testWithoutContact(floating_base_robot);
  testWithContact(floating_base_robot_contact);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}