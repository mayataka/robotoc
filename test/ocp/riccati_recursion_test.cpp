#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/hybrid/hybrid_container.hpp"

namespace idocp {

class RiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    N = 20;
    max_num_impulse = 5;
    nproc = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  using HybridKKTMatrix = hybrid_container<KKTMatrix, ImpulseKKTMatrix>;
  using HybridKKTResidual = hybrid_container<KKTResidual, ImpulseKKTResidual>;
  using HybridRiccatiFactorization = hybrid_container<RiccatiFactorization, RiccatiFactorization>;
  using HybridRiccatiRecursion = hybrid_container<RiccatiFactorizer, ImpulseRiccatiFactorizer>;

  HybridKKTMatrix createHybridKKTMatrix(const Robot& robot) const;
  HybridKKTResidual createHybridKKTResidual(const Robot& robot) const;

  DiscreteEvent createDiscreteEvent(const Robot& robot, const ContactStatus& pre_contact_status) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  static void RiccatiIsZero(const RiccatiFactorization& riccati);

  static void RiccatiIsSame(const RiccatiFactorization& lhs, const RiccatiFactorization& rhs);

  void testBackwardRiccatiRecursion(const Robot& robot) const;

  int N, max_num_impulse, nproc;
  double T, t, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


RiccatiRecursionTest::HybridKKTMatrix RiccatiRecursionTest::createHybridKKTMatrix(const Robot& robot) const {
  HybridKKTMatrix kkt_matrix = HybridKKTMatrix(N+1, KKTMatrix(robot), max_num_impulse, ImpulseKKTMatrix(robot));
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  for (int i=0; i<=N; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    tmp = Eigen::MatrixXd::Random(dimu, dimu);
    kkt_matrix[i].Quu() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimu, dimu);
    kkt_matrix[i].Qxu().setRandom();
    if (robot.has_floating_base()) {
      kkt_matrix[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix[i].Fvq().setRandom();
    kkt_matrix[i].Fvv().setRandom();
    kkt_matrix[i].Fvu().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.impulse[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    if (robot.has_floating_base()) {
      kkt_matrix.impulse[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.impulse[i].Fvq().setRandom();
    kkt_matrix.impulse[i].Fvv().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.aux[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    tmp = Eigen::MatrixXd::Random(dimu, dimu);
    kkt_matrix.aux[i].Quu() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimu, dimu);
    kkt_matrix.aux[i].Qxu().setRandom();
    if (robot.has_floating_base()) {
      kkt_matrix.aux[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.aux[i].Fvq().setRandom();
    kkt_matrix.aux[i].Fvv().setRandom();
    kkt_matrix.aux[i].Fvu().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.lift[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    tmp = Eigen::MatrixXd::Random(dimu, dimu);
    kkt_matrix.lift[i].Quu() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimu, dimu);
    kkt_matrix.lift[i].Qxu().setRandom();
    if (robot.has_floating_base()) {
      kkt_matrix.lift[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.lift[i].Fvq().setRandom();
    kkt_matrix.lift[i].Fvv().setRandom();
    kkt_matrix.lift[i].Fvu().setRandom();
  }
  return kkt_matrix;
}


RiccatiRecursionTest::HybridKKTResidual RiccatiRecursionTest::createHybridKKTResidual(const Robot& robot) const {
  HybridKKTResidual kkt_residual = HybridKKTResidual(N+1, KKTResidual(robot), max_num_impulse, ImpulseKKTResidual(robot));
  for (int i=0; i<=N; ++i) {
    kkt_residual[i].lx().setRandom();
    kkt_residual[i].lu().setRandom();
    kkt_residual[i].Fx().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    kkt_residual.impulse[i].lx().setRandom();
    kkt_residual.impulse[i].Fx().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    kkt_residual.aux[i].lx().setRandom();
    kkt_residual.aux[i].lu().setRandom();
    kkt_residual.aux[i].Fx().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    kkt_residual.lift[i].lx().setRandom();
    kkt_residual.lift[i].lu().setRandom();
    kkt_residual.lift[i].Fx().setRandom();
  }
  return kkt_residual;
}


DiscreteEvent RiccatiRecursionTest::createDiscreteEvent(const Robot& robot, const ContactStatus& pre_contact_status) const {
  DiscreteEvent discrete_event(robot);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  while (!discrete_event.existDiscreteEvent()) {
    post_contact_status.setRandom();
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


ContactSequence RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, T, N);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  for (int i=0; i<max_num_impulse; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    tmp.eventTime = i * 0.15 + 0.01 * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void RiccatiRecursionTest::RiccatiIsZero(const RiccatiFactorization& riccati) {
  EXPECT_TRUE(riccati.Pqq.isZero());
  EXPECT_TRUE(riccati.Pqv.isZero());
  EXPECT_TRUE(riccati.Pvq.isZero());
  EXPECT_TRUE(riccati.Pvv.isZero());
  EXPECT_TRUE(riccati.sq.isZero());
  EXPECT_TRUE(riccati.sv.isZero());
  EXPECT_TRUE(riccati.Pi.isIdentity()); // Default value of Pi is identity.
  EXPECT_TRUE(riccati.pi.isZero());
  EXPECT_TRUE(riccati.N.isZero());
  EXPECT_TRUE(riccati.n.isZero());
}


void RiccatiRecursionTest::RiccatiIsSame(const RiccatiFactorization& lhs, const RiccatiFactorization& rhs) {
  EXPECT_TRUE(lhs.Pqq.isApprox(rhs.Pqq));
  EXPECT_TRUE(lhs.Pqv.isApprox(rhs.Pqv));
  EXPECT_TRUE(lhs.Pvq.isApprox(rhs.Pvq));
  EXPECT_TRUE(lhs.Pvv.isApprox(rhs.Pvv));
  EXPECT_TRUE(lhs.sq.isApprox(rhs.sq));
  EXPECT_TRUE(lhs.sv.isApprox(rhs.sv));
  EXPECT_TRUE(lhs.Pi.isApprox(rhs.Pi));
  EXPECT_TRUE(lhs.pi.isApprox(rhs.pi));
  EXPECT_TRUE(lhs.N.isApprox(rhs.N));
  EXPECT_TRUE(lhs.n.isApprox(rhs.n));
}


void RiccatiRecursionTest::testBackwardRiccatiRecursion(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = createHybridKKTMatrix(robot);
  auto kkt_residual= createHybridKKTResidual(robot);
  HybridRiccatiFactorization factorization(N+1, RiccatiFactorization(robot), max_num_impulse, RiccatiFactorization(robot));
  RiccatiRecursion factorizer(robot, T, N, max_num_impulse, nproc);
  factorizer.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  for (int i=0; i<N; ++i) {
    RiccatiIsZero(factorization[i]);
  }
  EXPECT_TRUE(factorization[N].Pqq.isApprox(kkt_matrix[N].Qqq()));
  EXPECT_TRUE(factorization[N].Pqv.isZero());
  EXPECT_TRUE(factorization[N].Pvq.isZero());
  EXPECT_TRUE(factorization[N].Pvv.isApprox(kkt_matrix[N].Qvv()));
  EXPECT_TRUE(factorization[N].sq.isApprox(-1*kkt_residual[N].lq()));
  EXPECT_TRUE(factorization[N].sv.isApprox(-1*kkt_residual[N].lv()));
  for (int i=0; i<max_num_impulse; ++i) {
    RiccatiIsZero(factorization.impulse[i]);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    RiccatiIsZero(factorization.aux[i]);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    RiccatiIsZero(factorization.lift[i]);
  }
  factorizer.backwardRiccatiRecursion(contact_sequence, kkt_matrix, kkt_residual, factorization);

}


TEST_F(RiccatiRecursionTest, fixedBase) {
  testBackwardRiccatiRecursion(fixed_base_robot);
}


TEST_F(RiccatiRecursionTest, floating_base) {
  testBackwardRiccatiRecursion(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}