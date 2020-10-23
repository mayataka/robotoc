#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/impulse/impulse_riccati_matrix_factorizer.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"


namespace idocp {

class ImpulseRiccatiMatrixFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
  }

  virtual void TearDown() {
  }

  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(ImpulseRiccatiMatrixFactorizerTest, fixed_base) {
  const int dimv = fixed_base_robot_.dimv();
  const Eigen::MatrixXd P_seed = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  const Eigen::MatrixXd P_next = P_seed * P_seed.transpose();
  RiccatiFactorization riccati_next(fixed_base_robot_);
  riccati_next.Pqq = P_next.topLeftCorner(dimv, dimv);
  riccati_next.Pqv = P_next.topRightCorner(dimv, dimv);
  riccati_next.Pvq = P_next.bottomLeftCorner(dimv, dimv);
  riccati_next.Pvv = P_next.bottomRightCorner(dimv, dimv);
  const Eigen::VectorXd s_next = Eigen::VectorXd::Random(2*dimv);
  riccati_next.sq = s_next.head(dimv);
  riccati_next.sv = s_next.tail(dimv);
  ImpulseKKTMatrix kkt_matrix(fixed_base_robot_);
  const Eigen::MatrixXd Qxx_seed = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  kkt_matrix.Qxx() = Qxx_seed * Qxx_seed.transpose();
  kkt_matrix.Fqq = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  ImpulseKKTResidual kkt_residual(fixed_base_robot_);
  kkt_residual.Fq() = Eigen::VectorXd::Random(dimv);
  kkt_residual.Fv() = Eigen::VectorXd::Random(dimv);
  kkt_residual.lq() = Eigen::VectorXd::Random(dimv);
  kkt_residual.lv() = Eigen::VectorXd::Random(dimv);
  RiccatiFactorization riccati(fixed_base_robot_);
  ImpulseRiccatiMatrixFactorizer factorizer(fixed_base_robot_);
  factorizer.factorize(kkt_matrix, kkt_residual, riccati_next, riccati);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix.Fqq;
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq;
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv;
  const Eigen::MatrixXd P_ref = kkt_matrix.Qxx() + A.transpose() * P_next * A;
  const Eigen::VectorXd s_ref = A.transpose() * s_next 
                                - A.transpose() * P_next * kkt_residual.Fx() 
                                - kkt_residual.lx();
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}


TEST_F(ImpulseRiccatiMatrixFactorizerTest, floating_base) {
  const int dimv = floating_base_robot_.dimv();
  const Eigen::MatrixXd P_seed = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  const Eigen::MatrixXd P_next = P_seed * P_seed.transpose();
  RiccatiFactorization riccati_next(floating_base_robot_);
  riccati_next.Pqq = P_next.topLeftCorner(dimv, dimv);
  riccati_next.Pqv = P_next.topRightCorner(dimv, dimv);
  riccati_next.Pvq = P_next.bottomLeftCorner(dimv, dimv);
  riccati_next.Pvv = P_next.bottomRightCorner(dimv, dimv);
  const Eigen::VectorXd s_next = Eigen::VectorXd::Random(2*dimv);
  riccati_next.sq = s_next.head(dimv);
  riccati_next.sv = s_next.tail(dimv);
  ImpulseKKTMatrix kkt_matrix(floating_base_robot_);
  const Eigen::MatrixXd Qxx_seed = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  kkt_matrix.Qxx() = Qxx_seed * Qxx_seed.transpose();
  kkt_matrix.Fqq = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqq.topLeftCorner(6, 6).setRandom();
  kkt_matrix.Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  ImpulseKKTResidual kkt_residual(floating_base_robot_);
  kkt_residual.Fq() = Eigen::VectorXd::Random(dimv);
  kkt_residual.Fv() = Eigen::VectorXd::Random(dimv);
  kkt_residual.lq() = Eigen::VectorXd::Random(dimv);
  kkt_residual.lv() = Eigen::VectorXd::Random(dimv);
  RiccatiFactorization riccati(floating_base_robot_);
  ImpulseRiccatiMatrixFactorizer factorizer(floating_base_robot_);
  factorizer.factorize(kkt_matrix, kkt_residual, riccati_next, riccati);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix.Fqq;
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq;
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv;
  const Eigen::MatrixXd P_ref = kkt_matrix.Qxx() + A.transpose() * P_next * A;
  const Eigen::VectorXd s_ref = A.transpose() * s_next 
                                - A.transpose() * P_next * kkt_residual.Fx() 
                                - kkt_residual.lx();
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