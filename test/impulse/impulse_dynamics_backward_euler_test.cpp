#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/schur_complement.hpp"


namespace idocp {

class ImpulseDynamicsBackwardEulerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double t_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ImpulseDynamicsBackwardEulerTest, linearizeImpulseDynamicsBackwardEulerFixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.ldv = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.ldv = kkt_residual.ldv;
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  Eigen::MatrixXd dimd_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dimd_dq, dimd_ddv);
  robot.dRNEAPartialdFext(contact_status, dimd_df);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix_ref.Cq(), 
                                          kkt_matrix_ref.Cv());
  kkt_residual_ref.lq() += dimd_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv += dimd_ddv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dimd_df.transpose() * s.beta;
  EXPECT_TRUE(kkt_residual.dv_res.isApprox(kkt_residual_ref.dv_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.ldv.isApprox(kkt_residual_ref.ldv));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  const double l1norm_ref = kkt_residual_ref.dv_res.lpNorm<1>() 
                            + kkt_residual_ref.C().lpNorm<1>();
  const double squarednorm_ref = kkt_residual_ref.dv_res.squaredNorm()
                                  + kkt_residual_ref.C().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ImpulseDynamicsBackwardEulerTest, condenseImpulseDynamicsBackwardEulerFixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  const Eigen::VectorXd Qdvdv_diag = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::MatrixXd Qdvdv_ref = Qdvdv_diag.asDiagonal();
  const Eigen::VectorXd Qff_diag = Eigen::VectorXd::Random(contact_status.dimf()).array().abs();
  const Eigen::MatrixXd Qff_ref = Qff_diag.asDiagonal();
  kkt_matrix.Qdvdv = Qdvdv_ref;
  kkt_matrix.Qff() = Qff_ref;
  kkt_matrix_ref.Qdvdv = Qdvdv_ref;
  kkt_matrix_ref.Qff() = Qff_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.ldv = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.ldv = kkt_residual.ldv;
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.condenseImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dimd_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dimd_dq, dimd_ddv);
  robot.updateKinematics(s.q, s.v);
  robot.dRNEAPartialdFext(contact_status, dimd_df);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix_ref.Cq(), 
                                          kkt_matrix_ref.Cv());
  EXPECT_TRUE(kkt_matrix_ref.Cv().isApprox(-1*dimd_df.transpose()));
  kkt_residual_ref.lq() += dimd_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv += dimd_ddv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dimd_df.transpose() * s.beta;
  const Eigen::MatrixXd Minv 
      = dimd_ddv.llt().solve(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()));
  kkt_matrix_ref.Qqq() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_dq;
  kkt_matrix_ref.Qqf() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_df;
  kkt_residual_ref.lq() -= dimd_dq.transpose() * Minv * kkt_residual_ref.ldv;
  kkt_residual_ref.lq() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * kkt_residual_ref.dv_res;
  kkt_matrix_ref.Qff() += dimd_df.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_df;
  kkt_residual_ref.lf() -= dimd_df.transpose() * Minv * kkt_residual_ref.ldv;
  kkt_residual_ref.lf() += dimd_df.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * kkt_residual_ref.dv_res;
  kkt_matrix_ref.Fvq = - Minv * dimd_dq;
  kkt_matrix_ref.Fvf() = - Minv * dimd_df;
  kkt_residual_ref.Fv() = - Minv * kkt_residual.dv_res;
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(kkt_matrix_ref.Qff()));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(kkt_matrix_ref.Qfq()));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(kkt_matrix_ref.Qfv()));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(kkt_matrix_ref.Qqq()));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(kkt_matrix_ref.Qqv()));
  EXPECT_TRUE(kkt_matrix.Qvq().isApprox(kkt_matrix_ref.Qvq()));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(kkt_matrix_ref.Qvv()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Fvq.isApprox(kkt_matrix_ref.Fvq));
  EXPECT_TRUE(kkt_matrix.Fvv.isApprox(kkt_matrix_ref.Fvv));
  EXPECT_TRUE(kkt_matrix.Fvf().isApprox(kkt_matrix_ref.Fvf()));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, contact_status);
  id.computeCondensedDirection(kkt_matrix, kkt_residual, d);
  const Eigen::VectorXd ddv_ref = - Minv * dimd_dq * d.dq() 
                                  - Minv * dimd_df * d.df()
                                  - Minv * kkt_residual_ref.dv_res;
  EXPECT_TRUE(ddv_ref.isApprox(d.ddv));
  const Eigen::VectorXd dbeta_ref = Minv * kkt_matrix.Qdvdv * Minv * dimd_dq * d.dq() 
                                    + Minv * kkt_matrix.Qdvdv * Minv * dimd_df * d.df()
                                    - Minv * d.dgmm()
                                    - Minv * (kkt_residual_ref.ldv - kkt_matrix.Qdvdv * Minv * kkt_residual_ref.dv_res);
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
}


TEST_F(ImpulseDynamicsBackwardEulerTest, computeImpulseDynamicsResidualFixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.computeImpulseDynamicsResidual(robot, contact_status, s, kkt_residual);
  const double violation_ref = kkt_residual_ref.dv_res.lpNorm<1>() + kkt_residual_ref.C().lpNorm<1>();
  robot.updateKinematics(s.q, s.v);
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  const double l1norm_ref = kkt_residual_ref.dv_res.lpNorm<1>() 
                            + kkt_residual_ref.C().lpNorm<1>();
  const double squarednorm_ref = kkt_residual_ref.dv_res.squaredNorm()
                                  + kkt_residual_ref.C().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ImpulseDynamicsBackwardEulerTest, linearizeImpulseDynamicsBackwardEulerFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.ldv = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.ldv = kkt_residual.ldv;
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dimd_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dimd_dq, dimd_ddv);
  robot.dRNEAPartialdFext(contact_status, dimd_df);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix_ref.Cq(), 
                                          kkt_matrix_ref.Cv());
  kkt_residual_ref.lq() += dimd_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv += dimd_ddv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dimd_df.transpose() * s.beta;
  EXPECT_TRUE(kkt_residual.dv_res.isApprox(kkt_residual_ref.dv_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.ldv.isApprox(kkt_residual_ref.ldv));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  const double l1norm_ref = kkt_residual_ref.dv_res.lpNorm<1>() 
                            + kkt_residual_ref.C().lpNorm<1>();
  const double squarednorm_ref = kkt_residual_ref.dv_res.squaredNorm()
                                  + kkt_residual_ref.C().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ImpulseDynamicsBackwardEulerTest, condenseImpulseDynamicsBackwardEulerFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  const Eigen::VectorXd Qdvdv_diag = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::MatrixXd Qdvdv_ref = Qdvdv_diag.asDiagonal();
  const Eigen::VectorXd Qff_diag = Eigen::VectorXd::Random(contact_status.dimf()).array().abs();
  const Eigen::MatrixXd Qff_ref = Qff_diag.asDiagonal();
  kkt_matrix.Qdvdv = Qdvdv_ref;
  kkt_matrix.Qff() = Qff_ref;
  kkt_matrix_ref.Qdvdv = Qdvdv_ref;
  kkt_matrix_ref.Qff() = Qff_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.ldv = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.ldv = kkt_residual.ldv;
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.condenseImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dimd_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dimd_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dimd_dq, dimd_ddv);
  robot.updateKinematics(s.q, s.v);
  robot.dRNEAPartialdFext(contact_status, dimd_df);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix_ref.Cq(), 
                                          kkt_matrix_ref.Cv());
  EXPECT_TRUE(kkt_matrix_ref.Cv().isApprox(-1*dimd_df.transpose()));
  kkt_residual_ref.lq() += dimd_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv += dimd_ddv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dimd_df.transpose() * s.beta;
  const Eigen::MatrixXd Minv 
      = dimd_ddv.llt().solve(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()));
  kkt_matrix_ref.Qqq() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_dq;
  kkt_matrix_ref.Qqf() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_df;
  kkt_residual_ref.lq() -= dimd_dq.transpose() * Minv * kkt_residual_ref.ldv;
  kkt_residual_ref.lq() += dimd_dq.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * kkt_residual_ref.dv_res;
  kkt_matrix_ref.Qff() += dimd_df.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * dimd_df;
  kkt_residual_ref.lf() -= dimd_df.transpose() * Minv * kkt_residual_ref.ldv;
  kkt_residual_ref.lf() += dimd_df.transpose() * Minv * kkt_matrix_ref.Qdvdv * Minv * kkt_residual_ref.dv_res;
  kkt_matrix_ref.Fvq = - Minv * dimd_dq;
  kkt_matrix_ref.Fvf() = - Minv * dimd_df;
  kkt_residual_ref.Fv() = - Minv * kkt_residual.dv_res;
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(kkt_matrix_ref.Qff()));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(kkt_matrix_ref.Qfq()));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(kkt_matrix_ref.Qfv()));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(kkt_matrix_ref.Qqq()));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(kkt_matrix_ref.Qqv()));
  EXPECT_TRUE(kkt_matrix.Qvq().isApprox(kkt_matrix_ref.Qvq()));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(kkt_matrix_ref.Qvv()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Fvq.isApprox(kkt_matrix_ref.Fvq));
  EXPECT_TRUE(kkt_matrix.Fvv.isApprox(kkt_matrix_ref.Fvv));
  EXPECT_TRUE(kkt_matrix.Fvf().isApprox(kkt_matrix_ref.Fvf()));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, contact_status);
  id.computeCondensedDirection(kkt_matrix, kkt_residual, d);
  const Eigen::VectorXd ddv_ref = - Minv * dimd_dq * d.dq() 
                                  - Minv * dimd_df * d.df()
                                  - Minv * kkt_residual_ref.dv_res;
  EXPECT_TRUE(ddv_ref.isApprox(d.ddv));
  const Eigen::VectorXd dbeta_ref = Minv * kkt_matrix.Qdvdv * Minv * dimd_dq * d.dq() 
                                    + Minv * kkt_matrix.Qdvdv * Minv * dimd_df * d.df()
                                    - Minv * d.dgmm()
                                    - Minv * (kkt_residual_ref.ldv - kkt_matrix.Qdvdv * Minv * kkt_residual_ref.dv_res);
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
}


TEST_F(ImpulseDynamicsBackwardEulerTest, computeImpulseDynamicsBackwardEulerResidualFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  ImpulseKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  ImpulseKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.computeImpulseDynamicsResidual(robot, contact_status, s, kkt_residual);
  const double violation_ref = kkt_residual_ref.dv_res.lpNorm<1>() + kkt_residual_ref.C().lpNorm<1>();
  robot.updateKinematics(s.q, s.v);
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual_ref.dv_res);
  robot.computeContactVelocityResidual(contact_status, kkt_residual_ref.C());
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  const double l1norm_ref = kkt_residual_ref.dv_res.lpNorm<1>() 
                            + kkt_residual_ref.C().lpNorm<1>();
  const double squarednorm_ref = kkt_residual_ref.dv_res.squaredNorm()
                                  + kkt_residual_ref.C().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}