#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/unocp/backward_unriccati_recursion_factorizer.hpp"


namespace idocp {

class BackwardUnRiccatiRecursionFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    robot = Robot(urdf);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();

    unkkt_matrix = SplitUnKKTMatrix(robot);
    const Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
    unkkt_matrix.Q = seed * seed.transpose();

    unkkt_residual = SplitUnKKTResidual(robot);
    unkkt_residual.KKT_residual.setRandom();
  }

  virtual void TearDown() {
  }

  SplitRiccatiFactorization createRiccatiFactorization() const;

  std::string urdf;
  Robot robot;
  double dtau;
  int dimv;
  SplitUnKKTMatrix unkkt_matrix;
  SplitUnKKTResidual unkkt_residual;
};


SplitRiccatiFactorization BackwardUnRiccatiRecursionFactorizerTest::createRiccatiFactorization() const {
  SplitRiccatiFactorization riccati(robot, false);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  return riccati; 
}


TEST_F(BackwardUnRiccatiRecursionFactorizerTest, test) {
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization();
  const SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  const SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  BackwardUnRiccatiRecursionFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, dtau, unkkt_matrix, unkkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv).setIdentity();
  A.topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  A.bottomLeftCorner(dimv, dimv).setZero();
  A.bottomRightCorner(dimv, dimv).setIdentity();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimv);
  B.topRows(dimv).setZero();
  B.bottomRows(dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = unkkt_matrix_ref.Qxx() + A.transpose() * P_next * A;
  const Eigen::MatrixXd H_ref = unkkt_matrix_ref.Qax().transpose() + A.transpose() * P_next * B;
  const Eigen::MatrixXd G_ref = unkkt_matrix_ref.Qaa() + B.transpose() * P_next * B;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  const Eigen::VectorXd la_ref = B.transpose() * P_next * unkkt_residual_ref.Fx() - B.transpose() * sx_next + unkkt_residual_ref.la();
  EXPECT_TRUE(F_ref.isApprox(unkkt_matrix.Qxx()));
  EXPECT_TRUE(unkkt_matrix.Qxx().isApprox(unkkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(H_ref.isApprox(unkkt_matrix.Qax().transpose()));
  EXPECT_TRUE(G_ref.isApprox(unkkt_matrix.Qaa()));
  EXPECT_TRUE(unkkt_matrix.Qaa().isApprox(unkkt_matrix.Qaa().transpose()));
  EXPECT_TRUE(la_ref.isApprox(unkkt_residual.la()));
  SplitRiccatiFactorization riccati = createRiccatiFactorization();
  LQRStateFeedbackPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeRiccatiFactorization(riccati_next, unkkt_matrix, unkkt_residual, lqr_policy, dtau, riccati);
  const Eigen::MatrixXd P_ref = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * unkkt_residual_ref.Fx() - unkkt_residual_ref.lx() - H_ref * lqr_policy.k;
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}