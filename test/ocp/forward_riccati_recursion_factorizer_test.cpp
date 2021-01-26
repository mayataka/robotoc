#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/forward_riccati_recursion_factorizer.hpp"

namespace idocp {

class ForwardRiccatiRecursionFactorizerTest : public ::testing::Test {
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

  SplitKKTMatrix createKKTMatrix(const Robot& robot) const;
  SplitKKTResidual createKKTResidual(const Robot& robot) const;
  SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot) const;

  void test(const Robot& robot) const;

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


SplitKKTMatrix ForwardRiccatiRecursionFactorizerTest::createKKTMatrix(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq() = seed * seed.transpose();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Qvv() = seed * seed.transpose();
  kkt_matrix.Qqu().setRandom();
  kkt_matrix.Qvu().setRandom();
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  kkt_matrix.Quu() = seed * seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
    kkt_matrix.Fqv().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  return kkt_matrix;
}


SplitKKTResidual ForwardRiccatiRecursionFactorizerTest::createKKTResidual(const Robot& robot) const {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


SplitRiccatiFactorization ForwardRiccatiRecursionFactorizerTest::createRiccatiFactorization(const Robot& robot) const {
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
  const int dimx = 2 * robot.dimv();
  seed = Eigen::MatrixXd::Random(dimx, dimx);
  riccati.N = seed * seed.transpose();
  riccati.Pi.setRandom();
  riccati.pi.setRandom();
  return riccati; 
}


void ForwardRiccatiRecursionFactorizerTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const SplitKKTMatrix kkt_matrix = createKKTMatrix(robot);
  if (!robot.hasFloatingBase()) {
    ASSERT_TRUE(kkt_matrix.Fqq().isZero());
    ASSERT_TRUE(kkt_matrix.Fqv().isZero());
  }
  const SplitKKTResidual kkt_residual = createKKTResidual(robot);
  ForwardRiccatiRecursionFactorizer factorizer(robot);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv).setIdentity();
  if (robot.hasFloatingBase()) {
    A.topLeftCorner(robot.dim_passive(), robot.dim_passive()) 
        = kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive());
  }
  A.topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    A.topRightCorner(dimv, dimv).topLeftCorner(robot.dim_passive(), robot.dim_passive()) 
        = kkt_matrix.Fqv().topLeftCorner(robot.dim_passive(), robot.dim_passive());
  }
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimu);
  B.bottomRows(dimv) = kkt_matrix.Fvu();
  const SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_next(robot), riccati_next_ref(robot);
  factorizer.factorizeStateTransition(riccati, kkt_matrix, kkt_residual, dtau, riccati_next);
  riccati_next_ref.Pi = A * riccati.Pi;
  riccati_next_ref.pi = A * riccati.pi + kkt_residual.Fx();
  EXPECT_TRUE(riccati_next.Pi.isApprox(riccati_next_ref.Pi));
  EXPECT_TRUE(riccati_next.pi.isApprox(riccati_next_ref.pi));
  factorizer.factorizeStateConstraintFactorization(riccati, kkt_matrix, dtau, riccati_next);
  riccati_next_ref.N = A * riccati.N * A.transpose();
  EXPECT_TRUE(riccati_next.N.isApprox(riccati_next_ref.N));
  std::cout << riccati_next_ref.N - riccati_next.N << std::endl;
}


TEST_F(ForwardRiccatiRecursionFactorizerTest, fixedBase) {
  test(fixed_base_robot);
}


TEST_F(ForwardRiccatiRecursionFactorizerTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}