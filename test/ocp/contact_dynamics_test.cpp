#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

class ContactDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condenseContactDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix_ref.Qxx() = kkt_matrix.Qxx();
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix_ref.Qaaff() = kkt_matrix.Qaaff();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  std::cout << "aaa" << std::endl;
  cd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  Eigen::MatrixXd Minv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.computeMinv(dID_da, Minv);
  Eigen::MatrixXd Minv_dID_dqv = Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv());
  Minv_dID_dqv.leftCols(robot.dimv()) = Minv * dID_dq;
  Minv_dID_dqv.rightCols(robot.dimv()) = Minv * dID_dv;
  kkt_matrix_ref.Qxx() += Minv_dID_dqv.transpose() * kkt_matrix_ref.Qaaff() * Minv_dID_dqv;
  Eigen::MatrixXd IO = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimu());
  IO.topRows(robot.dimv()).setIdentity();
  kkt_matrix_ref.Qxu() -= Minv_dID_dqv.transpose() * kkt_matrix_ref.Qaaff() * Minv * IO;
  Eigen::VectorXd Minv_ID = Eigen::VectorXd::Zero(robot.dimv());
  Minv_ID = Minv * ID;
  kkt_residual_ref.lx() -= Minv_dID_dqv.transpose() * kkt_residual_ref.la;
  kkt_residual_ref.lx() += Minv_dID_dqv.transpose() * kkt_matrix_ref.Qaaff() * Minv_ID;
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix_ref.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxu().isApprox(kkt_matrix_ref.Qxu()));
  std::cout << kkt_matrix.Qxx() - kkt_matrix_ref.Qxx() << std::endl;
  std::cout << kkt_matrix.Qxu() - kkt_matrix_ref.Qxu() << std::endl;
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  cd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d.dgmm(), d);
  Eigen::VectorXd da_ref = - Minv * ID;
  da_ref -= Minv_dID_dqv * d.dx();
  da_ref += Minv * IO * d.du();
  EXPECT_TRUE(da_ref.isApprox(d.da()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


// TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFixedBaseWithoutContacts) {
//   std::vector<int> contact_frames = {18};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(fixed_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active = {false};
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ContactDynamics rd(robot);
//   rd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
// }


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


// TEST_F(ContactDynamicsTest, condenseContactDynamicsFixedBaseWithContact) {
//   std::vector<int> contact_frames = {18};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(fixed_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active = {true};
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
//   kkt_matrix.Quu = Quu_ref;
//   kkt_matrix_ref.Quu = Quu_ref;
//   kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
//   kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual_ref.lq() = kkt_residual.lq();
//   kkt_residual_ref.lv() = kkt_residual.lv();
//   kkt_residual_ref.la() = kkt_residual.la();
//   kkt_residual_ref.lf() = kkt_residual.lf();
//   kkt_residual_ref.lu = kkt_residual.lu;
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
//   Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.dRNEAPartialdFext(contact_status, du_df);
//   robot.computeBaumgarteResidual(contact_status, dtau_, dtau_, kkt_residual_ref.C());
//   robot.computeBaumgarteDerivatives(contact_status, dtau_, dtau_, kkt_matrix_ref.Cq(), 
//                                     kkt_matrix_ref.Cv(), kkt_matrix_ref.Ca());
//   kkt_residual_ref.lu -= dtau_ * s.beta;
//   Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
//   kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
//                           + du_dq.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
//   kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
//                           + du_dv.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
//   kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
//                           + du_da.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
//   kkt_residual_ref.lf() = dtau_ * du_df.transpose() * s.beta 
//                           + du_df.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
//   EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
//   EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
//   EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
//   EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
//   EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
//   EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
//   EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
//   EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(contact_status.dimf())));
//   EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_df.leftCols(contact_status.dimf())));
//   EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
//   EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
//   EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
//   EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
//   EXPECT_TRUE(kkt_matrix.Cf().isApprox(kkt_matrix_ref.Cf()));
//   SplitDirection d = SplitDirection::Random(robot, contact_status);
//   rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
//   Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
//   du_ref += du_dq * d.dq();
//   du_ref += du_dv * d.dv();
//   du_ref += du_da * d.da();
//   du_ref += du_df.leftCols(contact_status.dimf()) * d.df();
//   EXPECT_TRUE(du_ref.isApprox(d.du));
//   Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
//   EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
//   std::cout << "du_dq" << std::endl;
//   std::cout << du_dq << std::endl;
//   std::cout << "du_dv" << std::endl;
//   std::cout << du_dv << std::endl;
//   std::cout << "du_da" << std::endl;
//   std::cout << du_da << std::endl;
//   std::cout << "du_df" << std::endl;
//   std::cout << du_df << std::endl;
//   std::cout << "Cq" << std::endl;
//   std::cout << kkt_matrix.Cq() << std::endl;
//   std::cout << "Cv" << std::endl;
//   std::cout << kkt_matrix.Cv() << std::endl;
//   std::cout << "Ca" << std::endl;
//   std::cout << kkt_matrix.Ca() << std::endl;
//   std::cout << "Cf" << std::endl;
//   std::cout << kkt_matrix.Cf() << std::endl;
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + kkt_residual_ref.C().head(contact_status.dimf()).lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                                   + kkt_residual_ref.C().head(contact_status.dimf()).squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
//   const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
//   EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq+du_df*df_dq));
//   EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv+du_df*df_dv));
// }


// TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFixedBaseWithContacts) {
//   std::vector<int> contact_frames = {18};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(fixed_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active = {true};
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   robot.computeBaumgarteResidual(contact_status, dtau_, dtau_, 
//                                  kkt_residual_ref.C_contacts());
//   kkt_residual_ref.u_res -= s.u;
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + kkt_residual_ref.C().head(contact_status.dimf()).lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                                   + kkt_residual_ref.C().head(contact_status.dimf()).squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
// }


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(false);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(robot.dim_passive()) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(kkt_residual_ref.lu_passive));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>() + dtau_ * s.u_passive.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm() + dtau_ * dtau_ * s.u_passive.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


// TEST_F(ContactDynamicsTest, condenseContactDynamicsFloatingBaseWithoutContacts) {
//   std::vector<int> contact_frames = {14, 24, 34, 44};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(floating_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active;
//   for (const auto frame : contact_frames) {
//     is_contact_active.push_back(false);
//   }
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
//   kkt_matrix.Quu = Quu_ref;
//   kkt_matrix_ref.Quu = Quu_ref;
//   kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual_ref.lq() = kkt_residual.lq();
//   kkt_residual_ref.lv() = kkt_residual.lv();
//   kkt_residual_ref.la() = kkt_residual.la();
//   kkt_residual_ref.lu = kkt_residual.lu;
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
//   Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
//   Eigen::MatrixXd Cu 
//       = Eigen::MatrixXd::Zero(robot.dim_passive()+contact_status.dimf(), robot.dimv());  
//   Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
//   kkt_residual_ref.C().head(robot.dim_passive()) 
//       = dtau_ * s.u.head(robot.dim_passive());
//   kkt_residual_ref.lu -= dtau_ * s.beta;
//   kkt_residual_ref.lu += Cu.transpose() * s.mu_stack();
//   Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
//   kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
//                           + du_dq.transpose() * lu_condensed
//                           + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
//   kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
//                           + du_dv.transpose() * lu_condensed
//                           + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
//   kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
//                           + du_da.transpose() * lu_condensed
//                           + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
//   kkt_residual_ref.C() += Cu * kkt_residual_ref.u_res;
//   kkt_matrix_ref.Cq() += Cu * du_dq;
//   kkt_matrix_ref.Cv() += Cu * du_dv;
//   kkt_matrix_ref.Ca() += Cu * du_da;
//   kkt_matrix_ref.Cf() += Cu * du_df;
//   EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
//   EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
//   EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
//   EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
//   EXPECT_TRUE(kkt_residual.lf().isZero());
//   EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
//   EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
//   EXPECT_TRUE(kkt_matrix.Qaf().isZero());
//   EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qff().isZero());
//   EXPECT_TRUE(kkt_matrix.Qfq().isZero());
//   EXPECT_TRUE(kkt_matrix.Qfv().isZero());
//   EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
//   EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
//   EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
//   EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
//   EXPECT_TRUE(kkt_matrix.Cf().isZero());
//   SplitDirection d = SplitDirection::Random(robot);
//   rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
//   Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
//   du_ref += du_dq * d.dq();
//   du_ref += du_dv * d.dv();
//   du_ref += du_da * d.da();
//   EXPECT_TRUE(du_ref.isApprox(d.du));
//   Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
//   dbeta_ref += Cu.transpose() * d.dmu() / dtau_;
//   EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
//   std::cout << "du_dq" << std::endl;
//   std::cout << du_dq << std::endl;
//   std::cout << "du_dv" << std::endl;
//   std::cout << du_dv << std::endl;
//   std::cout << "du_da" << std::endl;
//   std::cout << du_da << std::endl;
//   std::cout << "du_df" << std::endl;
//   std::cout << du_df << std::endl;
//   std::cout << "Cq" << std::endl;
//   std::cout << kkt_matrix.Cq() << std::endl;
//   std::cout << "Cv" << std::endl;
//   std::cout << kkt_matrix.Cv() << std::endl;
//   std::cout << "Ca" << std::endl;
//   std::cout << kkt_matrix.Ca() << std::endl;
//   std::cout << "Cf" << std::endl;
//   std::cout << kkt_matrix.Cf() << std::endl;
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + dtau_ * s.u.head(6).lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                             + dtau_ * dtau_ * s.u.head(6).squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
//   const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
//   EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq));
//   EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv));
// }


// TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFloatingBaseWithoutContacts) {
//   std::vector<int> contact_frames = {14, 24, 34, 44};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(floating_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active;
//   for (const auto frame : contact_frames) {
//     is_contact_active.push_back(false);
//   }
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + dtau_ * s.u.head(6).lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                             + dtau_ * dtau_ * s.u.head(6).squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
// }


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(robot.dim_passive()) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(kkt_residual_ref.lu_passive));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + s.u_passive.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + s.u_passive.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


// TEST_F(ContactDynamicsTest, condenseContactDynamicsFloatingBaseWithContacts) {
//   std::vector<int> contact_frames = {14, 24, 34, 44};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(floating_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active;
//   for (const auto frame : contact_frames) {
//     is_contact_active.push_back(rnd()%2==0);
//   }
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
//   kkt_matrix.Quu = Quu_ref;
//   kkt_matrix_ref.Quu = Quu_ref;
//   kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
//   kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
//   kkt_residual_ref.lq() = kkt_residual.lq();
//   kkt_residual_ref.lv() = kkt_residual.lv();
//   kkt_residual_ref.la() = kkt_residual.la();
//   kkt_residual_ref.lf() = kkt_residual.lf();
//   kkt_residual_ref.lu = kkt_residual.lu;
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
//   Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());  
//   Eigen::MatrixXd Cu 
//       = Eigen::MatrixXd::Zero(robot.dim_passive()+contact_status.dimf(), robot.dimv());  
//   Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
//   robot.setContactForces(contact_status, s.f);
//   robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
//   kkt_residual_ref.u_res -= s.u;
//   robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
//   robot.updateKinematics(s.q, s.v, s.a);
//   robot.dRNEAPartialdFext(contact_status, du_df);
//   robot.computeBaumgarteResidual(contact_status, dtau_, dtau_, kkt_residual_ref.C().tail(contact_status.dimf()));
//   kkt_residual_ref.C().head(robot.dim_passive()) 
//       = dtau_ * s.u.head(robot.dim_passive());
//   robot.computeBaumgarteDerivatives(contact_status, dtau_, dtau_, kkt_matrix_ref.Cq().bottomRows(contact_status.dimf()), 
//                                     kkt_matrix_ref.Cv().bottomRows(contact_status.dimf()), 
//                                     kkt_matrix_ref.Ca().bottomRows(contact_status.dimf()));
//   kkt_residual_ref.lu -= dtau_ * s.beta;
//   kkt_residual_ref.lu += Cu.transpose() * s.mu_stack();
//   Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
//   kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
//                           + du_dq.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
//   kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
//                           + du_dv.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
//   kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
//                           + du_da.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
//   kkt_residual_ref.lf() = dtau_ * du_df.transpose() * s.beta 
//                           + du_df.transpose() * lu_condensed 
//                           + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
//   kkt_residual_ref.C() += Cu * kkt_residual_ref.u_res;
//   kkt_matrix_ref.Cq() += Cu * du_dq;
//   kkt_matrix_ref.Cv() += Cu * du_dv;
//   kkt_matrix_ref.Ca() += Cu * du_da;
//   kkt_matrix_ref.Cf() += Cu * du_df;
//   EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
//   EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
//   EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
//   EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
//   EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
//   EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
//   EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
//   EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(contact_status.dimf())));
//   EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_df.leftCols(contact_status.dimf())));
//   EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(contact_status.dimf()).transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
//   EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
//   EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
//   EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
//   EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
//   EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
//   EXPECT_TRUE(kkt_matrix.Cf().isApprox(kkt_matrix_ref.Cf()));
//   SplitDirection d = SplitDirection::Random(robot, contact_status);
//   rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
//   Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
//   du_ref += du_dq * d.dq();
//   du_ref += du_dv * d.dv();
//   du_ref += du_da * d.da();
//   if (contact_status.hasActiveContacts()) {
//     du_ref += du_df.leftCols(contact_status.dimf()) * d.df();
//   }
//   EXPECT_TRUE(du_ref.isApprox(d.du));
//   Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
//   dbeta_ref += Cu.transpose() * d.dmu() / dtau_;
//   EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
//   std::cout << "du_dq" << std::endl;
//   std::cout << du_dq << std::endl;
//   std::cout << "du_dv" << std::endl;
//   std::cout << du_dv << std::endl;
//   std::cout << "du_da" << std::endl;
//   std::cout << du_da << std::endl;
//   std::cout << "du_df" << std::endl;
//   std::cout << du_df << std::endl;
//   std::cout << "Cq" << std::endl;
//   std::cout << kkt_matrix.Cq() << std::endl;
//   std::cout << "Cv" << std::endl;
//   std::cout << kkt_matrix.Cv() << std::endl;
//   std::cout << "Ca" << std::endl;
//   std::cout << kkt_matrix.Ca() << std::endl;
//   std::cout << "Cf" << std::endl;
//   std::cout << kkt_matrix.Cf() << std::endl;
//   std::cout << "Cu" << std::endl;
//   std::cout << Cu << std::endl;
//   std::cout << "Cu.transpose()" << std::endl;
//   std::cout << Cu.transpose() << std::endl;
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
//   const double l1norm_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>() 
//                             + dtau_ * s.u.head(6).lpNorm<1>()
//                             + kkt_residual_ref.C_contacts().lpNorm<1>();
//   const double squarednorm_ref = dtau_ * dtau_ * kkt_residual_ref.u_res.squaredNorm()
//                             + dtau_ * dtau_ * s.u.head(6).squaredNorm()
//                             + kkt_residual_ref.C_contacts().squaredNorm();
//   EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
//   EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
//   const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
//   const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(contact_status.dimf(), robot.dimv());
//   Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
//   rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
//   EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq+du_df*df_dq));
//   EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv+du_df*df_dv));
// }


// TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFloatingBaseWithContacts) {
//   std::vector<int> contact_frames = {14, 24, 34, 44};
//   ContactStatus contact_status(contact_frames.size());
//   Robot robot(floating_base_urdf_, contact_frames);
//   std::random_device rnd;
//   std::vector<bool> is_contact_active;
//   for (const auto frame : contact_frames) {
//     is_contact_active.push_back(rnd()%2==0);
//   }
//   contact_status.setContactStatus(is_contact_active);
//   SplitSolution s = SplitSolution::Random(robot, contact_status);
//   KKTResidual kkt_residual(robot);
//   kkt_residual.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix(robot);
//   kkt_matrix.setContactStatus(contact_status);
//   KKTResidual kkt_residual_ref(robot);
//   kkt_residual_ref.setContactStatus(contact_status);
//   KKTMatrix kkt_matrix_ref(robot);
//   kkt_matrix_ref.setContactStatus(contact_status);
//   ContactDynamics rd(robot);
//   robot.updateKinematics(s.q, s.v, s.a);
//   rd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
//   const double l1norm = rd.l1NormContactDynamicsResidual(dtau_, kkt_residual);
//   const double squarednorm = rd.squaredNormContactDynamicsResidual(dtau_, kkt_residual);
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