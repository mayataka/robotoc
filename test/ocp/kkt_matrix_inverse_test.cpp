#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/kkt_matrix_inverse.hpp"


namespace idocp {

class KKTMatrixInverseTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(KKTMatrixInverseTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.setContactStatus(robot);
  KKTMatrixInverse matrix(robot);
  matrix.setContactStatus(robot);
  EXPECT_EQ(matrix.dimKKT(), composition.dimKKT());
  EXPECT_EQ(matrix.max_dimKKT(), composition.max_dimKKT());
  const int dimKKT = matrix.dimKKT();
  Eigen::MatrixXd matrix_ref = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  matrix.KKT_matrix_inverse() = matrix_ref;
  EXPECT_TRUE(matrix.KKT_matrix_inverse().isApprox(matrix_ref));
  const int dimx = 2*robot.dimv();
  EXPECT_TRUE(matrix.auxiliaryMatrix().isApprox(matrix_ref.topLeftCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.backwardCorrectionSerialCoeff().isApprox(matrix_ref.topRightCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.backwardCorrectionParallelCoeff().isApprox(matrix_ref.bottomRightCorner(dimKKT-dimx, dimx)));
  EXPECT_TRUE(matrix.forwardCorrectionSerialCoeff().isApprox(matrix_ref.bottomLeftCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.forwardCorrectionParallelCoeff().isApprox(matrix_ref.topLeftCorner(dimKKT-dimx, dimx)));
  matrix.setZero();
  EXPECT_TRUE(matrix.KKT_matrix_inverse().isZero());
}


TEST_F(KKTMatrixInverseTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.setContactStatus(robot);
  KKTMatrixInverse matrix(robot);
  matrix.setContactStatus(robot);
  EXPECT_EQ(matrix.dimKKT(), composition.dimKKT());
  EXPECT_EQ(matrix.max_dimKKT(), composition.max_dimKKT());
  const int dimKKT = matrix.dimKKT();
  Eigen::MatrixXd matrix_ref = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  matrix.KKT_matrix_inverse() = matrix_ref;
  EXPECT_TRUE(matrix.KKT_matrix_inverse().isApprox(matrix_ref));
  const int dimx = 2*robot.dimv();
  EXPECT_TRUE(matrix.auxiliaryMatrix().isApprox(matrix_ref.topLeftCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.backwardCorrectionSerialCoeff().isApprox(matrix_ref.topRightCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.backwardCorrectionParallelCoeff().isApprox(matrix_ref.bottomRightCorner(dimKKT-dimx, dimx)));
  EXPECT_TRUE(matrix.forwardCorrectionSerialCoeff().isApprox(matrix_ref.bottomLeftCorner(dimx, dimx)));
  EXPECT_TRUE(matrix.forwardCorrectionParallelCoeff().isApprox(matrix_ref.topLeftCorner(dimKKT-dimx, dimx)));
  matrix.setZero();
  EXPECT_TRUE(matrix.KKT_matrix_inverse().isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}