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
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/schur_complement.hpp"


namespace idocp {

class ImpulseDynamicsForwardEulerTest : public ::testing::Test {
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


TEST_F(ImpulseDynamicsForwardEulerTest, linearizeImpulseDynamicsForwardEulerFixedBase) {
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
  ImpulseDynamicsForwardEuler id(robot);
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


TEST_F(ImpulseDynamicsForwardEulerTest, condenseImpulseDynamicsForwardEulerFixedBase) {
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
  ImpulseDynamicsForwardEuler id(robot);
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
  SchurComplement schur_complement(robot.dimv(), contact_status.dimf());
  Eigen::MatrixXd MJTJinv = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                  robot.dimv()+contact_status.dimf());
  schur_complement.invertWithZeroBottomRightCorner(robot.dimv(), contact_status.dimf(),
                                                   dimd_ddv, kkt_matrix_ref.Cv(),
                                                   MJTJinv);
  Eigen::MatrixXd dimdc_dqv = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                    2*robot.dimv());
  dimdc_dqv.topLeftCorner(robot.dimv(), robot.dimv()) = dimd_dq;
  dimdc_dqv.bottomLeftCorner(contact_status.dimf(), robot.dimv()) = kkt_matrix_ref.Cq();
  dimdc_dqv.bottomRightCorner(contact_status.dimf(), robot.dimv()) = kkt_matrix_ref.Cv();
  const Eigen::MatrixXd MJTJinv_dimdc_dqv = MJTJinv * dimdc_dqv;
  Eigen::MatrixXd Qdvdvff = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                  robot.dimv()+contact_status.dimf());
  Qdvdvff.topLeftCorner(robot.dimv(), robot.dimv()) = kkt_matrix_ref.Qdvdv;
  Qdvdvff.bottomRightCorner(contact_status.dimf(), contact_status.dimf()) = kkt_matrix_ref.Qff();
  kkt_matrix_ref.Qxx() += MJTJinv_dimdc_dqv.transpose() * Qdvdvff.diagonal().asDiagonal() * MJTJinv_dimdc_dqv;
  Eigen::VectorXd ldvf = Eigen::VectorXd::Zero(robot.dimv()+contact_status.dimf());
  ldvf.head(robot.dimv()) = kkt_residual_ref.ldv;
  ldvf.tail(contact_status.dimf()) = - kkt_residual_ref.lf();
  Eigen::VectorXd imdc = Eigen::VectorXd::Zero(robot.dimv()+contact_status.dimf());
  imdc.head(robot.dimv()) = kkt_residual_ref.dv_res;
  imdc.tail(contact_status.dimf()) = kkt_residual_ref.C();
  kkt_residual_ref.lx() -= MJTJinv_dimdc_dqv.transpose() * ldvf;
  kkt_residual_ref.lx() += MJTJinv_dimdc_dqv.transpose() * Qdvdvff.diagonal().asDiagonal() * MJTJinv * imdc;
  kkt_matrix_ref.Fvq = - MJTJinv_dimdc_dqv.topLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - MJTJinv_dimdc_dqv.topRightCorner(robot.dimv(), robot.dimv());
  kkt_residual_ref.Fv() = - (MJTJinv * imdc).head(robot.dimv());
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(kkt_matrix_ref.Qqq()));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(kkt_matrix_ref.Qqv()));
  EXPECT_TRUE(kkt_matrix.Qvq().isApprox(kkt_matrix_ref.Qvq()));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(kkt_matrix_ref.Qvv()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Fvq.isApprox(kkt_matrix_ref.Fvq));
  EXPECT_TRUE(kkt_matrix.Fvv.isApprox(kkt_matrix_ref.Fvv));
  std::cout << kkt_residual.lq() - kkt_residual_ref.lq() << std::endl;
  std::cout << kkt_residual.lv() - kkt_residual_ref.lv() << std::endl;
  // SplitDirection d = SplitDirection::Random(robot, contact_status);
  // rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  // Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
  // du_ref += du_dq * d.dq();
  // du_ref += du_dv * d.dv();
  // du_ref += du_da * d.da();
  // du_ref += du_df.leftCols(contact_status.dimf()) * d.df();
  // EXPECT_TRUE(du_ref.isApprox(d.du));
  // Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
  // EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  // std::cout << "du_dq" << std::endl;
  // std::cout << du_dq << std::endl;
  // std::cout << "du_dv" << std::endl;
  // std::cout << du_dv << std::endl;
  // std::cout << "du_da" << std::endl;
  // std::cout << du_da << std::endl;
  // std::cout << "du_df" << std::endl;
  // std::cout << du_df << std::endl;
  // std::cout << "Cq" << std::endl;
  // std::cout << kkt_matrix.Cq() << std::endl;
  // std::cout << "Cv" << std::endl;
  // std::cout << kkt_matrix.Cv() << std::endl;
  // std::cout << "Ca" << std::endl;
  // std::cout << kkt_matrix.Ca() << std::endl;
  // std::cout << "Cf" << std::endl;
  // std::cout << kkt_matrix.Cf() << std::endl;
  // const double l1norm = rd.l1NormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
  // const double squarednorm = rd.squaredNormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
  // const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
  //                           + kkt_residual_ref.C().head(contact_status.dimf()).lpNorm<1>();
  // const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
  //                                 + kkt_residual_ref.C().head(contact_status.dimf()).squaredNorm();
  // EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  // EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
  // const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  // const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  // const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
  // const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
  // Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  // rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
  // EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq+du_df*df_dq));
  // EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv+du_df*df_dv));
}


// TEST_F(ImpulseDynamicsForwardEulerTest, computeImpulseDynamicsForwardEulerResidualFixedBase) {
//   std::vector<int> contact_frames = {18};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(fixed_base_urdf_, contact_frames);
//   std::vector<bool> is_contact_active = {true};
//   contact_status.setContactStatus(is_contact_active);
//   ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
//   ImpulseKKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   ImpulseKKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   ImpulseKKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   ImpulseKKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ImpulseDynamicsForwardEuler rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.computeImpulseDynamicsForwardEulerResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   robot.computeBaumgarteResidual(contact_status, dtau_, dtau_, 
//                                  kkt_residual_ref.C_contacts());
//   kkt_residual_ref.u_res -= s.u;
//   const double l1norm = rd.l1NormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + kkt_residual_ref.C().head(contact_status.dimf()).lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                                   + kkt_residual_ref.C().head(contact_status.dimf()).squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
// }


TEST_F(ImpulseDynamicsForwardEulerTest, linearizeImpulseDynamicsForwardEulerFloatingBaseWithContacts) {
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
  ImpulseDynamicsForwardEuler id(robot);
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


TEST_F(ImpulseDynamicsForwardEulerTest, condenseImpulseDynamicsForwardEulerFloatingBaseWithContacts) {
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
  ImpulseDynamicsForwardEuler id(robot);
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
  SchurComplement schur_complement(robot.dimv(), robot.max_dimf());
  Eigen::MatrixXd MJTJinv = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                  robot.dimv()+contact_status.dimf());
  schur_complement.invertWithZeroBottomRightCorner(robot.dimv(), contact_status.dimf(),
                                                   dimd_ddv, kkt_matrix_ref.Cv(),
                                                   MJTJinv);
  Eigen::MatrixXd dimdc_dqv = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                    2*robot.dimv());
  dimdc_dqv.topLeftCorner(robot.dimv(), robot.dimv()) = dimd_dq;
  dimdc_dqv.bottomLeftCorner(contact_status.dimf(), robot.dimv()) = kkt_matrix_ref.Cq();
  dimdc_dqv.bottomRightCorner(contact_status.dimf(), robot.dimv()) = kkt_matrix_ref.Cv();
  const Eigen::MatrixXd MJTJinv_dimdc_dqv = MJTJinv * dimdc_dqv;
  Eigen::MatrixXd Qdvdvff = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), 
                                                  robot.dimv()+contact_status.dimf());
  Qdvdvff.topLeftCorner(robot.dimv(), robot.dimv()) = kkt_matrix_ref.Qdvdv;
  Qdvdvff.bottomRightCorner(contact_status.dimf(), contact_status.dimf()) = kkt_matrix_ref.Qff();
  kkt_matrix_ref.Qxx() += MJTJinv_dimdc_dqv.transpose() * Qdvdvff.diagonal().asDiagonal() * MJTJinv_dimdc_dqv;
  Eigen::VectorXd ldvf = Eigen::VectorXd::Zero(robot.dimv()+contact_status.dimf());
  ldvf.head(robot.dimv()) = kkt_residual_ref.ldv;
  ldvf.tail(contact_status.dimf()) = - kkt_residual_ref.lf();
  Eigen::VectorXd imdc = Eigen::VectorXd::Zero(robot.dimv()+contact_status.dimf());
  imdc.head(robot.dimv()) = kkt_residual_ref.dv_res;
  imdc.tail(contact_status.dimf()) = kkt_residual_ref.C();
  kkt_residual_ref.lx() -= MJTJinv_dimdc_dqv.transpose() * ldvf;
  kkt_residual_ref.lx() += MJTJinv_dimdc_dqv.transpose() * Qdvdvff.diagonal().asDiagonal() * MJTJinv * imdc;
  kkt_matrix_ref.Fvq = - MJTJinv_dimdc_dqv.topLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - MJTJinv_dimdc_dqv.topRightCorner(robot.dimv(), robot.dimv());
  kkt_residual_ref.Fv() = - (MJTJinv * imdc).head(robot.dimv());
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(kkt_matrix_ref.Qqq()));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(kkt_matrix_ref.Qqv()));
  EXPECT_TRUE(kkt_matrix.Qvq().isApprox(kkt_matrix_ref.Qvq()));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(kkt_matrix_ref.Qvv()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Fvq.isApprox(kkt_matrix_ref.Fvq));
  EXPECT_TRUE(kkt_matrix.Fvv.isApprox(kkt_matrix_ref.Fvv));
}


// TEST_F(ImpulseDynamicsForwardEulerTest, computeImpulseDynamicsForwardEulerResidualFloatingBaseWithContacts) {
//   std::vector<int> contact_frames = {14, 24, 34, 44};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(floating_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active;
//   for (const auto frame : contact_frames) {
//     is_contact_active.push_back(rnd()%2==0);
//   }
//   contact_status.setContactStatus(is_contact_active);
//   ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, contact_status);
//   ImpulseKKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   ImpulseKKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   ImpulseKKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   ImpulseKKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ImpulseDynamicsForwardEuler rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.computeImpulseDynamicsForwardEulerResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double l1norm = rd.l1NormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormImpulseDynamicsForwardEulerResidual(dtau_, kkt_residual);
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   robot.computeBaumgarteResidual(contact_status, dtau_, dtau_, kkt_residual_ref.C_contacts());
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + dtau_ * s.u.head(6).lpNorm<1>()
//                             + kkt_residual_ref.C_contacts().lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                             + dtau_ * dtau_ * s.u.head(6).squaredNorm()
//                             + kkt_residual_ref.C_contacts().squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
// }

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}