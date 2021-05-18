#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseBackwardRiccatiRecursionFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  ImpulseSplitKKTMatrix createKKTMatrix(const Robot& robot) const;
  ImpulseSplitKKTResidual createKKTResidual(const Robot& robot) const;
  SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot) const;

  void test(const Robot& robot) const;
};


ImpulseSplitKKTMatrix ImpulseBackwardRiccatiRecursionFactorizerTest::createKKTMatrix(const Robot& robot) const {
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq() = seed * seed.transpose();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Qvv() = seed * seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  return kkt_matrix;
}


ImpulseSplitKKTResidual ImpulseBackwardRiccatiRecursionFactorizerTest::createKKTResidual(const Robot& robot) const {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.lx.setRandom();
  kkt_residual.Fx.setRandom();
  return kkt_residual;
}


SplitRiccatiFactorization ImpulseBackwardRiccatiRecursionFactorizerTest::createRiccatiFactorization(const Robot& robot) const {
  SplitRiccatiFactorization riccati(robot);
  const int dimv = robot.dimv();
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


void ImpulseBackwardRiccatiRecursionFactorizerTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const auto riccati_next = createRiccatiFactorization(robot);
  auto kkt_matrix = createKKTMatrix(robot);
  auto kkt_residual = createKKTResidual(robot);
  const auto kkt_matrix_ref = kkt_matrix;
  const auto kkt_residual_ref = kkt_residual;
  ImpulseBackwardRiccatiRecursionFactorizer factorizer(robot);
  if (!robot.hasFloatingBase()) {
    ASSERT_TRUE(kkt_matrix.Fqq().isZero());
  }
  ASSERT_TRUE(kkt_matrix.Fqv().isZero());
  factorizer.factorizeKKTMatrix(riccati_next, kkt_matrix);
  if (!robot.hasFloatingBase()) {
    ASSERT_TRUE(kkt_matrix.Fqq().isZero());
  }
  ASSERT_TRUE(kkt_matrix.Fqv().isZero());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  if (robot.hasFloatingBase()) {
    A.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Fqq();
  }
  else {
    A.topLeftCorner(dimv, dimv).setIdentity();
  }
  A.topRightCorner(dimv, dimv).setZero();
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix_ref.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix_ref.Fvv();
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * P_next * A;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  SplitRiccatiFactorization riccati(robot), riccati_ref(robot);
  if (!robot.hasFloatingBase()) {
    ASSERT_TRUE(kkt_matrix.Fqq().isZero());
  }
  ASSERT_TRUE(kkt_matrix.Fqv().isZero());
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, riccati);
  if (!robot.hasFloatingBase()) {
    ASSERT_TRUE(kkt_matrix.Fqq().isZero());
  }
  ASSERT_TRUE(kkt_matrix.Fqv().isZero());
  const Eigen::MatrixXd P_ref = F_ref;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * kkt_residual_ref.Fx - kkt_residual_ref.lx;
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}


TEST_F(ImpulseBackwardRiccatiRecursionFactorizerTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(ImpulseBackwardRiccatiRecursionFactorizerTest, floating_base) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}