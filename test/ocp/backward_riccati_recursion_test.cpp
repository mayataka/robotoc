#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion.hpp"


namespace idocp {

class BackwardRiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


void BackwardRiccatiRecursionTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiFactorization riccati_next(robot);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pqq = seed * seed.transpose();
  riccati_next.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvv = seed * seed.transpose();
  riccati_next.sq.setRandom();
  riccati_next.sv.setRandom();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = seed * seed.transpose();
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Qqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = seed * seed.transpose();
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd Qqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Qvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Quu = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Qqu() = Qqu;
  kkt_matrix.Qvu() = Qvu;
  kkt_matrix.Quu() = Quu;
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  const KKTMatrix kkt_matrix_ref = kkt_matrix;
  const KKTResidual kkt_residual_ref = kkt_residual;
  BackwardRiccatiRecursion factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix, kkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Fqq();
  A.topRightCorner(dimv, dimv) = kkt_matrix_ref.Fqv();
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix_ref.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix_ref.Fvv();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimu);
  B.bottomRows(dimv) = kkt_matrix_ref.Fvu();
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx() + A.transpose() * P_next * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu() + A.transpose() * P_next * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Quu() + B.transpose() * P_next * B;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  const Eigen::VectorXd lu_ref = B.transpose() * P_next * kkt_residual_ref.Fx() - B.transpose() * sx_next + kkt_residual_ref.lu();
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
  EXPECT_TRUE(lu_ref.isApprox(kkt_residual.lu()));
  RiccatiFactorization riccati(robot), riccati_ref(robot);
  LQRStateFeedbackPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, dtau, riccati);
  const Eigen::MatrixXd P_ref = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * kkt_residual_ref.Fx() - kkt_residual_ref.lx() - H_ref * lqr_policy.k;
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}


TEST_F(BackwardRiccatiRecursionTest, fixedBase) {
  test(fixed_base_robot);
}


TEST_F(BackwardRiccatiRecursionTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}