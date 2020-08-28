#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"


namespace idocp {

class RiccatiMatrixFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(RiccatiMatrixFactorizerTest, fixed_base) {
  const int dimv = fixed_base_robot_.dimv();
  RiccatiMatrixFactorizer factorizer(fixed_base_robot_);
  Eigen::MatrixXd Pqq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pvq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pvv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq_ref = Qqq + Pqq;
  const Eigen::MatrixXd Qqv_ref = Qqv + dtau_ * Pqq + Pqv;
  const Eigen::MatrixXd Qvq_ref = Qvq + dtau_ * Pqq + Pvq;
  const Eigen::MatrixXd Qvv_ref = Qvv + dtau_ * dtau_ * Pqq 
                                      + dtau_ * (Pqv + Pvq) + Pvv;
  factorizer.factorizeF(dtau_, Pqq, Pqv, Pvq, Pvv, Qqq, Qqv, Qvq, Qvv);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  std::cout << "Qqq error:" << std::endl;
  std::cout << Qqq - Qqq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qqv error:" << std::endl;
  std::cout << Qqv - Qqv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qvq error:" << std::endl;
  std::cout << Qvq - Qvq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qvv error:" << std::endl;
  std::cout << Qvv - Qvv_ref << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqa_ref = Qqa + dtau_ * Pqv;
  const Eigen::MatrixXd Qva_ref = Qva + dtau_ * dtau_ * Pqv + dtau_ * Pvv;
  factorizer.factorizeH(dtau_, Pqv, Pvv, Qqa, Qva);
  EXPECT_TRUE(Qqa.isApprox(Qqa_ref));
  EXPECT_TRUE(Qva.isApprox(Qva_ref));
  std::cout << "Qqa error:" << std::endl;
  std::cout << Qqa - Qqa_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qva error:" << std::endl;
  std::cout << Qva - Qva_ref << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qaa_ref = Qaa + dtau_ * dtau_ * Pvv;
  factorizer.factorizeG(dtau_, Pvv, Qaa);
  EXPECT_TRUE(Qaa.isApprox(Qaa_ref));
  std::cout << "Qaa error:" << std::endl;
  std::cout << Qaa - Qaa_ref << std::endl;
  std::cout << std::endl;

  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd sv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd la_ref = la + dtau_ * Pvq * Fq + dtau_ * Pvv * Fv - dtau_ * sv;
  factorizer.factorize_la(dtau_, Pvq, Pvv, Fq, Fv, sv, la);

  const Eigen::MatrixXd Pqq_ref = Pqq;
  const Eigen::MatrixXd Pqv_ref = Pqv;
  const Eigen::MatrixXd Pvq_ref = Pvq;
  const Eigen::MatrixXd Pvv_ref = Pvv;
  Eigen::VectorXd sq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd sq_ref = sq;
  factorizer.correctP(Pqq, Pqv);
  EXPECT_TRUE(Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(Pqv.isApprox(Pqv_ref));
  factorizer.correct_s(sq);
  EXPECT_TRUE(sq.isApprox(sq_ref));
}


TEST_F(RiccatiMatrixFactorizerTest, floating_base) {
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  RiccatiMatrixFactorizer factorizer(floating_base_robot_);
  Eigen::MatrixXd Pqq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pvq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Pvv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Zero(dimq);
  floating_base_robot_.generateFeasibleConfiguration(q_prev);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(dimq);
  floating_base_robot_.generateFeasibleConfiguration(q);
  Eigen::VectorXd q_next = Eigen::VectorXd::Zero(dimq);
  floating_base_robot_.generateFeasibleConfiguration(q_next);
  Eigen::MatrixXd Fqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Fqq_prev = Eigen::MatrixXd::Zero(dimv, dimv);
  floating_base_robot_.dSubtractdConfigurationPlus(q, q_next, Fqq);
  floating_base_robot_.dSubtractdConfigurationMinus(q_prev, q, Fqq_prev);
  const Eigen::MatrixXd Qqq_ref = Qqq + Fqq.transpose() * Pqq * Fqq;
  const Eigen::MatrixXd Qqv_ref = Qqv + Fqq.transpose() * (dtau_ * Pqq + Pqv);
  const Eigen::MatrixXd Qvq_ref = Qvq + (dtau_ * Pqq + Pvq) * Fqq;
  const Eigen::MatrixXd Qvv_ref = Qvv + (dtau_*dtau_) * Pqq 
                                      + dtau_ * (Pvq + Pqv) + Pvv;
  factorizer.setStateEquationDerivatives(Fqq, Fqq_prev);
  factorizer.factorizeF(dtau_, Pqq, Pqv, Pvq, Pvv, Qqq, Qqv, Qvq, Qvv);
  EXPECT_TRUE(Qqq.isApprox(Qqq_ref));
  EXPECT_TRUE(Qqv.isApprox(Qqv_ref));
  EXPECT_TRUE(Qvq.isApprox(Qvq_ref));
  EXPECT_TRUE(Qvv.isApprox(Qvv_ref));
  std::cout << "Qqq error:" << std::endl;
  std::cout << Qqq - Qqq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qqv error:" << std::endl;
  std::cout << Qqv - Qqv_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qvq error:" << std::endl;
  std::cout << Qvq - Qvq_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qvv error:" << std::endl;
  std::cout << Qvv - Qvv_ref << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd Qqa = Eigen::MatrixXd::Random(dimv, dimv);
  Eigen::MatrixXd Qva = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqa_ref = Qqa + dtau_ * Fqq.transpose() * Pqv;
  const Eigen::MatrixXd Qva_ref = Qva + dtau_ * dtau_ * Pqv + dtau_ * Pvv;
  factorizer.factorizeH(dtau_, Pqv, Pvv, Qqa, Qva);
  EXPECT_TRUE(Qqa.isApprox(Qqa_ref));
  EXPECT_TRUE(Qva.isApprox(Qva_ref));
  std::cout << "Qqa error:" << std::endl;
  std::cout << Qqa - Qqa_ref << std::endl;
  std::cout << std::endl;
  std::cout << "Qva error:" << std::endl;
  std::cout << Qva - Qva_ref << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qaa_ref = Qaa + dtau_ * dtau_ * Pvv;
  factorizer.factorizeG(dtau_, Pvv, Qaa);
  EXPECT_TRUE(Qaa.isApprox(Qaa_ref));
  std::cout << "Qaa error:" << std::endl;
  std::cout << Qaa - Qaa_ref << std::endl;
  std::cout << std::endl;

  Eigen::VectorXd la = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd sv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd la_ref = la + dtau_ * Pvq * Fq + dtau_ * Pvv * Fv - dtau_ * sv;
  factorizer.factorize_la(dtau_, Pvq, Pvv, Fq, Fv, sv, la);

  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P.topLeftCorner(dimv, dimv) = Pqq;
  P.topRightCorner(dimv, dimv) = Pqv;
  P.bottomLeftCorner(dimv, dimv) = Pvq;
  P.bottomRightCorner(dimv, dimv) = Pvv;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = - Fqq_prev.transpose();
  const Eigen::MatrixXd A_trans_inv = A.inverse();
  std::cout << "A_trans_inv" << std::endl;
  std::cout << A_trans_inv.topLeftCorner(dimv, dimv) << std::endl;
  const Eigen::MatrixXd P_ref = A_trans_inv * P;
  Eigen::VectorXd sq = Eigen::VectorXd::Random(dimv);
  Eigen::VectorXd s = Eigen::VectorXd::Zero(2*dimv);
  s.head(dimv) = sq;
  s.tail(dimv) = sv;
  const Eigen::VectorXd s_ref = A_trans_inv * s;
  factorizer.correctP(Pqq, Pqv);
  EXPECT_TRUE(Pqq.isApprox(P_ref.topLeftCorner(dimv, dimv)));
  EXPECT_TRUE(Pqv.isApprox(P_ref.topRightCorner(dimv, dimv)));
  EXPECT_TRUE(Pvq.isApprox(P_ref.bottomLeftCorner(dimv, dimv)));
  EXPECT_TRUE(Pvv.isApprox(P_ref.bottomRightCorner(dimv, dimv)));
  factorizer.correct_s(sq);
  EXPECT_TRUE(sq.isApprox(s_ref.head(dimv)));
  EXPECT_TRUE(sv.isApprox(s_ref.tail(dimv)));

  std::cout << "Pqq - P_ref.topLeftCorner(dimv, dimv)" << std::endl;
  std::cout << Pqq - P_ref.topLeftCorner(dimv, dimv) << std::endl;
  std::cout << "Pqv - P_ref.topRightCorner(dimv, dimv)" << std::endl;
  std::cout << Pqv - P_ref.topRightCorner(dimv, dimv) << std::endl;
  std::cout << "sq - s_ref.head(dimv)" << std::endl;
  std::cout << sq - s_ref.head(dimv) << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}