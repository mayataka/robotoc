#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/switching_constraint_jacobian.hpp"
#include "robotoc/ocp/switching_constraint_residual.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/backward_riccati_recursion_factorizer.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace robotoc {

class RiccatiFactorizerTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  virtual void TearDown() {
  }
  double dt;
};


TEST_P(RiccatiFactorizerTest, backwardRecursion) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  const bool sto = true;
  bool has_next_sto_phase = true;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, lqr_policy, sto, has_next_sto_phase);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  backward_recursion_ref.factorizeHamiltonian(riccati_next, kkt_matrix_ref, riccati_ref, has_next_sto_phase);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.llt().solve(Eigen::MatrixXd::Identity(dimu, dimu));
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  lqr_policy_ref.T = - Ginv  * riccati_ref.psi_u;
  lqr_policy_ref.W = - Ginv  * riccati_ref.phi_u;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  backward_recursion_ref.factorizeSTOFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, 
                                                   lqr_policy_ref, riccati_ref, has_next_sto_phase);
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(lqr_policy_ref.isApprox(lqr_policy));
  has_next_sto_phase = false;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, lqr_policy, sto, has_next_sto_phase);
  EXPECT_TRUE(lqr_policy.W.isZero());
}


TEST_P(RiccatiFactorizerTest, backwardRecursionPhaseTransition) {
  const auto robot = GetParam();
  const double max_dts0 = std::abs(Eigen::VectorXd::Random(1)[0]);
  RiccatiFactorizer factorizer(robot, max_dts0);
  STOPolicy sto_policy(robot), sto_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  auto riccati_m = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_m_ref = riccati_m;
  bool sto = false;
  factorizer.backwardRiccatiRecursionPhaseTransition(riccati, riccati_m, sto_policy, sto);
  riccati_m_ref.P = riccati.P;
  riccati_m_ref.s = riccati.s;
  riccati_m_ref.Psi.setZero();
  riccati_m_ref.Phi = riccati.Psi;
  riccati_m_ref.xi = 0.0;
  riccati_m_ref.chi = 0.0;
  riccati_m_ref.rho = riccati.xi;
  riccati_m_ref.eta = 0.0;
  riccati_m_ref.iota = riccati.eta;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  sto = true;
  factorizer.backwardRiccatiRecursionPhaseTransition(riccati, riccati_m, sto_policy, sto);
  const double keps = std::sqrt(std::numeric_limits<double>::epsilon());
  double xi = riccati.xi - 2.0 * riccati.chi + riccati.rho;
  if (xi*max_dts0 < std::abs(riccati.eta-riccati.iota) || xi < keps) {
    xi = std::abs(xi) + std::abs(riccati.eta-riccati.iota) / max_dts0;
  }
  // STO policy
  sto_policy_ref.dtsdx  = - (1.0/xi) * (riccati.Psi-riccati.Phi);
  sto_policy_ref.dtsdts =   (1.0/xi) * (riccati.xi-riccati.chi);
  sto_policy_ref.dts0   = - (1.0/xi) * (riccati.eta-riccati.iota);
  riccati_m_ref.s.noalias()   
      += (1.0/xi) * (riccati.Psi-riccati.Phi) * (riccati.eta-riccati.iota);
  riccati_m_ref.Phi.noalias() 
      -= (1.0/xi) * (riccati.Psi-riccati.Phi) * (riccati.xi-riccati.chi);
  riccati_m_ref.rho
      = riccati.xi - (1.0/xi) * (riccati.xi-riccati.chi) * (riccati.xi-riccati.chi);
  riccati_m_ref.iota
      = riccati.eta - (1.0/xi) * (riccati.xi-riccati.chi)  * (riccati.eta-riccati.iota);
  EXPECT_TRUE(sto_policy.isApprox(sto_policy_ref));
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati_m.isApprox(riccati_m_ref));
}


TEST_P(RiccatiFactorizerTest, backwardRecursionWithSwitchingConstraint) {
  const auto robot = GetParam();
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimi = impulse_status.dimi();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  SwitchingConstraintJacobian sc_jacobian(robot);
  SwitchingConstraintResidual sc_residual(robot);
  sc_jacobian.setImpulseStatus(impulse_status);
  sc_residual.setImpulseStatus(impulse_status);
  sc_jacobian.Phix().setRandom();
  sc_jacobian.Phia().setRandom();
  sc_jacobian.Phiu().setRandom();
  sc_jacobian.Phit().setRandom();
  sc_residual.P().setRandom();

  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  SplitConstrainedRiccatiFactorization c_riccati(robot);
  bool sto = true;
  bool has_next_sto_phase = true;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, 
                                      sc_jacobian, sc_residual, 
                                      riccati, c_riccati, lqr_policy, sto, has_next_sto_phase);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);

  Eigen::MatrixXd GDtD = Eigen::MatrixXd::Zero(dimu+dimi, dimu+dimi);
  GDtD.topLeftCorner(dimu, dimu) = kkt_matrix.Quu;
  GDtD.topRightCorner(dimu, dimi) = sc_jacobian.Phiu().transpose();
  GDtD.bottomLeftCorner(dimi, dimu) = sc_jacobian.Phiu();
  const Eigen::MatrixXd GDtDinv = GDtD.inverse();
  Eigen::MatrixXd HtC = Eigen::MatrixXd::Zero(dimu+dimi, 2*dimv);
  HtC.topRows(dimu) = kkt_matrix_ref.Qxu.transpose();
  HtC.bottomRows(dimi) = sc_jacobian.Phix();
  Eigen::VectorXd htc = Eigen::VectorXd::Zero(dimu+dimi);
  htc.head(dimu) = kkt_residual_ref.lu;
  htc.tail(dimi) = sc_residual.P();
  const Eigen::MatrixXd KM = - GDtDinv * HtC;
  const Eigen::VectorXd km = - GDtDinv * htc;
  lqr_policy_ref.K = KM.topRows(dimu);
  lqr_policy_ref.k = km.head(dimu);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  Eigen::MatrixXd ODtD = GDtD;
  ODtD.topLeftCorner(dimu, dimu).setZero();
  const Eigen::MatrixXd MtKtODtDKM = KM.transpose() * ODtD * KM;
  riccati_ref.P -= MtKtODtDKM;
  riccati_ref.s -= sc_jacobian.Phix().transpose() * km.tail(dimi);
  backward_recursion_ref.factorizeHamiltonian(riccati_next, kkt_matrix_ref, riccati_ref,
                                              has_next_sto_phase);
  Eigen::VectorXd psict = Eigen::VectorXd::Zero(dimu+dimi);
  psict.head(dimu) = riccati_ref.psi_u;
  psict.tail(dimi) = sc_jacobian.Phit();
  const Eigen::VectorXd Tmt = - GDtDinv * psict;
  Eigen::VectorXd phict = Eigen::VectorXd::Zero(dimu+dimi);
  phict.head(dimu) = riccati_ref.phi_u;
  const Eigen::VectorXd Wmt_next = - GDtDinv * phict;
  lqr_policy_ref.T = Tmt.head(dimu);
  lqr_policy_ref.W = Wmt_next.head(dimu);
  const Eigen::MatrixXd M_ref = KM.bottomRows(dimi);
  const Eigen::VectorXd m_ref = km.tail(dimi);
  const Eigen::VectorXd mt_ref = Tmt.tail(dimi);
  const Eigen::VectorXd mt_next_ref = Wmt_next.tail(dimi);

  EXPECT_TRUE(c_riccati.M().isApprox(M_ref));
  EXPECT_TRUE(c_riccati.m().isApprox(m_ref));
  EXPECT_TRUE(c_riccati.mt().isApprox(mt_ref));
  EXPECT_TRUE(c_riccati.mt_next().isApprox(mt_next_ref));

  backward_recursion_ref.factorizeSTOFactorization(riccati_next, kkt_matrix, 
                                                   kkt_residual, lqr_policy_ref, 
                                                   riccati_ref, has_next_sto_phase);
  riccati_ref.Psi += M_ref.transpose() * sc_jacobian.Phit();
  riccati_ref.xi += mt_ref.dot(sc_jacobian.Phit());
  if (has_next_sto_phase) {
    riccati_ref.chi += mt_next_ref.dot(sc_jacobian.Phit());
  }
  else {
    riccati_ref.chi = 0.0;
  }
  riccati_ref.eta += m_ref.dot(sc_jacobian.Phit());

  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(lqr_policy_ref.K.isApprox(lqr_policy.K));
  EXPECT_TRUE(lqr_policy_ref.k.isApprox(lqr_policy.k));
  if (!lqr_policy_ref.T.isZero()) {
    EXPECT_TRUE(lqr_policy_ref.T.isApprox(lqr_policy.T)); // if near zero, approx does not work well
  }
  if (!lqr_policy_ref.W.isZero()) {
    EXPECT_TRUE(lqr_policy_ref.W.isApprox(lqr_policy.W)); // if near zero, approx does not work well
  }
  std::cout << "impulse_status: " << impulse_status << std::endl;

  // Tests when has_next_sto_phase = false
  has_next_sto_phase = false;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, 
                                      sc_jacobian, sc_residual, 
                                      riccati, c_riccati, lqr_policy, sto, has_next_sto_phase);
  EXPECT_TRUE(lqr_policy.W.isZero());
  EXPECT_TRUE(c_riccati.mt_next().isZero());
}


TEST_P(RiccatiFactorizerTest, backwardRecursionImpulse) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateImpulseSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateImpulseSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  bool sto = true;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, sto);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
  backward_recursion_ref.factorizeSTOFactorization(riccati_next, kkt_matrix_ref, 
                                                   kkt_residual_ref, riccati_ref);
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
}


TEST_P(RiccatiFactorizerTest, computeSwitchingTimeDirection) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  auto sto_policy = STOPolicy(robot);
  sto_policy.dtsdx.setRandom();
  sto_policy.dtsdts = Eigen::VectorXd::Random(1)[0];
  sto_policy.dts0 = Eigen::VectorXd::Random(1)[0];
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  bool has_next_sto_phase = false;
  RiccatiFactorizer::computeSwitchingTimeDirection(sto_policy, d, has_next_sto_phase);
  d_ref.dts_next = sto_policy.dtsdx.dot(d.dx) + sto_policy.dts0;
  EXPECT_TRUE(d.isApprox(d_ref));
  has_next_sto_phase = true;
  RiccatiFactorizer::computeSwitchingTimeDirection(sto_policy, d, has_next_sto_phase);
  d_ref.dts_next += sto_policy.dtsdts * d.dts;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_P(RiccatiFactorizerTest, forwardRecursion) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  const Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.llt().solve(Eigen::MatrixXd::Identity(dimu, dimu));
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  auto d = SplitDirection::Random(robot);
  auto d_next = SplitDirection::Random(robot);
  auto d_ref = d;
  auto d_next_ref = d_next;
  bool sto = false;
  bool sto_next = false;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, lqr_policy, 
                                     d, d_next, sto, sto_next);
  d_ref.du = lqr_policy_ref.K * d_ref.dx + lqr_policy_ref.k;
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx;
  d_next_ref.dv() += kkt_matrix_ref.Fvu * d_ref.du;
  d_next_ref.dts = d_ref.dts;
  d_next_ref.dts_next = d_ref.dts_next;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  sto = true;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, lqr_policy, 
                                     d, d_next, sto, sto_next);
  d_ref.du += lqr_policy_ref.T * (d_ref.dts_next-d_ref.dts);
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx + kkt_matrix_ref.fx * (d_ref.dts_next-d_ref.dts);
  d_next_ref.dv() += kkt_matrix_ref.Fvu * d_ref.du;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  sto_next = true;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, lqr_policy, 
                                     d, d_next, sto, sto_next);
  d_ref.du -= lqr_policy_ref.W * d_ref.dts_next;
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx + kkt_matrix_ref.fx * (d_ref.dts_next-d_ref.dts);
  d_next_ref.dv() += kkt_matrix_ref.Fvu * d_ref.du;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  sto = false;
  sto_next = false;
  factorizer.computeCostateDirection(riccati, d, sto, sto_next);
  d_ref.dlmdgmm = riccati.P * d.dx - riccati.s;
  EXPECT_TRUE(d.isApprox(d_ref));
  sto = true;
  factorizer.computeCostateDirection(riccati, d, sto, sto_next);
  d_ref.dlmdgmm += riccati.Psi * (d_ref.dts_next-d_ref.dts);
  EXPECT_TRUE(d.isApprox(d_ref));
  sto_next = true;
  factorizer.computeCostateDirection(riccati, d, sto, sto_next);
  d_ref.dlmdgmm -= riccati.Phi * d_ref.dts_next;
  EXPECT_TRUE(d.isApprox(d_ref));
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  SplitConstrainedRiccatiFactorization c_riccati(robot);
  c_riccati.setImpulseStatus(impulse_status.dimi());
  c_riccati.M().setRandom();
  c_riccati.m().setRandom();
  d.setImpulseStatus(impulse_status);
  d.dxi().setRandom();
  d_ref = d;
  sto = false;
  sto_next = false;
  factorizer.computeLagrangeMultiplierDirection(c_riccati, d, sto, sto_next);
  d_ref.dxi() = c_riccati.M() * d.dx + c_riccati.m();
  EXPECT_TRUE(d.isApprox(d_ref));
  sto = true;
  factorizer.computeLagrangeMultiplierDirection(c_riccati, d, sto, sto_next);
  d_ref.dxi() += c_riccati.mt() * (d_ref.dts_next-d_ref.dts);
  EXPECT_TRUE(d.isApprox(d_ref));
  sto_next = true;
  factorizer.computeLagrangeMultiplierDirection(c_riccati, d, sto, sto_next);
  d_ref.dxi() -= c_riccati.mt_next() * d_ref.dts_next;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_P(RiccatiFactorizerTest, forwardRecursionImpulse) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix = testhelper::CreateImpulseSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateImpulseSplitKKTResidual(robot);
  RiccatiFactorizer factorizer(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  auto d = ImpulseSplitDirection::Random(robot);
  auto d_next = SplitDirection::Random(robot);
  auto d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d, d_next);
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx;
  d_next_ref.dts = d.dts;
  d_next_ref.dts_next = d.dts_next;
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  auto d_ref = d;
  bool sto = false;
  factorizer.computeCostateDirection(riccati, d, sto);
  d_ref.dlmdgmm = riccati.P * d.dx - riccati.s;
  EXPECT_TRUE(d.isApprox(d_ref));
  sto = true;
  factorizer.computeCostateDirection(riccati, d, sto);
  d_ref.dlmdgmm -= riccati.Phi * d_ref.dts_next;
  EXPECT_TRUE(d.isApprox(d_ref));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, RiccatiFactorizerTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateQuadrupedalRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}