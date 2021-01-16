#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

class SplitImpulseKKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void testSize(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testInverse(const Robot& robot, const ImpulseStatus& impulse_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitImpulseKKTMatrixTest::testSize(const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTMatrix matrix(robot);
  matrix.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = impulse_status.dimf();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimKKT(), 4*dimv+2*dimf);
  EXPECT_EQ(matrix.Fqf().rows(), dimv);
  EXPECT_EQ(matrix.Fqf().cols(), dimf);
  EXPECT_EQ(matrix.Fqq().rows(), dimv);
  EXPECT_EQ(matrix.Fqq().cols(), dimv);
  EXPECT_EQ(matrix.Fqv().rows(), dimv);
  EXPECT_EQ(matrix.Fqv().cols(), dimv);
  EXPECT_EQ(matrix.Fvf().rows(), dimv);
  EXPECT_EQ(matrix.Fvf().cols(), dimf);
  EXPECT_EQ(matrix.Fvq().rows(), dimv);
  EXPECT_EQ(matrix.Fvq().cols(), dimv);
  EXPECT_EQ(matrix.Fvv().rows(), dimv);
  EXPECT_EQ(matrix.Fvv().cols(), dimv);
  EXPECT_EQ(matrix.Fxf().rows(), 2*dimv);
  EXPECT_EQ(matrix.Fxf().cols(), dimf);
  EXPECT_EQ(matrix.Fxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Fxx().cols(), 2*dimv);
  EXPECT_EQ(matrix.Pq().rows(), dimf);
  EXPECT_EQ(matrix.Pq().cols(), dimv);
  EXPECT_EQ(matrix.Vq().rows(), dimf);
  EXPECT_EQ(matrix.Vq().cols(), dimv);
  EXPECT_EQ(matrix.Vv().rows(), dimf);
  EXPECT_EQ(matrix.Vv().cols(), dimv);
  EXPECT_EQ(matrix.Qdvdvff().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qdvdvff().cols(), dimv+dimf);
  EXPECT_EQ(matrix.Qdvdv().rows(), dimv);
  EXPECT_EQ(matrix.Qdvdv().cols(), dimv);
  EXPECT_EQ(matrix.Qff().rows(), dimf);
  EXPECT_EQ(matrix.Qff().cols(), dimf);
  EXPECT_EQ(matrix.Qfq().rows(), dimf);
  EXPECT_EQ(matrix.Qfq().cols(), dimv);
  EXPECT_EQ(matrix.Qqq().rows(), dimv);
  EXPECT_EQ(matrix.Qqq().cols(), dimv);
  EXPECT_EQ(matrix.Qqv().rows(), dimv);
  EXPECT_EQ(matrix.Qqv().cols(), dimv);
  EXPECT_EQ(matrix.Qvq().cols(), dimv);
  EXPECT_EQ(matrix.Qvv().rows(), dimv);
  EXPECT_EQ(matrix.Qvv().cols(), dimv);
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_EQ(matrix.Fqq_prev.rows(), dimv);
  EXPECT_EQ(matrix.Fqq_prev.cols(), dimv);
  const Eigen::MatrixXd Fqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Vq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Vv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qdvdv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Fqf() = Fqf;
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvf() = Fvf;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
  matrix.Pq() = Pq;
  matrix.Vq() = Vq;
  matrix.Vv() = Vv;
  matrix.Qdvdv() = Qdvdv;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  EXPECT_TRUE(matrix.Fqf().isApprox(Fqf));
  EXPECT_TRUE(matrix.Fqq().isApprox(Fqq));
  EXPECT_TRUE(matrix.Fqv().isApprox(Fqv));
  EXPECT_TRUE(matrix.Fvf().isApprox(Fvf));
  EXPECT_TRUE(matrix.Fvq().isApprox(Fvq));
  EXPECT_TRUE(matrix.Fvv().isApprox(Fvv));
  EXPECT_TRUE(matrix.Fxf().topRows(dimv).isApprox(Fqf));
  EXPECT_TRUE(matrix.Fxf().bottomRows(dimv).isApprox(Fvf));
  EXPECT_TRUE(matrix.Fxx().topLeftCorner(dimv, dimv).isApprox(Fqq));
  EXPECT_TRUE(matrix.Fxx().topRightCorner(dimv, dimv).isApprox(Fqv));
  EXPECT_TRUE(matrix.Fxx().bottomLeftCorner(dimv, dimv).isApprox(Fvq));
  EXPECT_TRUE(matrix.Fxx().bottomRightCorner(dimv, dimv).isApprox(Fvv));
  EXPECT_TRUE(matrix.Pq().isApprox(Pq));
  EXPECT_TRUE(matrix.Vq().isApprox(Vq));
  EXPECT_TRUE(matrix.Vv().isApprox(Vv));
  EXPECT_TRUE(matrix.Qdvdv().isApprox(Qdvdv));
  EXPECT_TRUE(matrix.Qff().isApprox(Qff));
  EXPECT_TRUE(matrix.Qdvdvff().topLeftCorner(dimv, dimv).isApprox(Qdvdv));
  EXPECT_TRUE(matrix.Qdvdvff().bottomRightCorner(dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.Qfq().isApprox(Qfq));
  EXPECT_TRUE(matrix.Qqq().isApprox(Qqq));
  EXPECT_TRUE(matrix.Qqv().isApprox(Qqv));
  EXPECT_TRUE(matrix.Qvq().isApprox(Qvq));
  EXPECT_TRUE(matrix.Qvv().isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxx().topLeftCorner(dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().topRightCorner(dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().bottomLeftCorner(dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().bottomRightCorner(dimv, dimv).isApprox(Qvv));
}


void SplitImpulseKKTMatrixTest::testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTMatrix matrix(robot);
  matrix.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  const Eigen::MatrixXd Fqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Pq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Vq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Vv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qdvdv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Fqf() = Fqf;
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvf() = Fvf;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
  matrix.Pq() = Pq;
  matrix.Vq() = Vq;
  matrix.Vv() = Vv;
  matrix.Qdvdv() = Qdvdv;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  ImpulseSplitKKTMatrix matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Fxx().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  if (impulse_status.hasActiveImpulse()) {
    matrix_ref.Fxf().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
    matrix_ref.Pq().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
    matrix_ref.Vq().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
    matrix_ref.Vv().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
  }
  else {
    matrix_ref.Fxf().setRandom();
    matrix_ref.Pq().setRandom();
    matrix_ref.Vq().setRandom();
    matrix_ref.Vv().setRandom();
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
  }
  matrix_ref.Qdvdvff().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  if (dimf > 0) {
    matrix_ref.Qfq().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
    matrix_ref.Qqf().setRandom();
    EXPECT_FALSE(matrix.isApprox(matrix_ref));
    matrix_ref = matrix;
    EXPECT_TRUE(matrix.isApprox(matrix_ref));
  }
  matrix_ref.Qxx().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
}


TEST_F(SplitImpulseKKTMatrixTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  ImpulseStatus impulse_status(contact_frames.size());
  impulse_status.setImpulseStatus({false});
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(SplitImpulseKKTMatrixTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_impulse_active = {false, false, false, false};
  ImpulseStatus impulse_status(contact_frames.size());
  impulse_status.setImpulseStatus(is_impulse_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  is_impulse_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_impulse_active.push_back(rnd()%2==0);
  }
  impulse_status.setImpulseStatus(is_impulse_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}