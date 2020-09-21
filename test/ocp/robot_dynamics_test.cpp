#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/robot_dynamics.hpp"


namespace idocp {

class RobotDynamicsTest : public ::testing::Test {
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


TEST_F(RobotDynamicsTest, augmentRobotDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {false};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  rd.augmentRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual_ref.lq() += dtau_ * du_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * du_dv.transpose() * s.beta; 
  kkt_residual_ref.la() += dtau_ * du_da.transpose() * s.beta;
  kkt_residual_ref.lu -= dtau_ * s.beta;
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_residual.C().isZero());
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
}


TEST_F(RobotDynamicsTest, condenseRobotDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {false};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
  kkt_matrix.Quu = Quu_ref;
  kkt_matrix_ref.Quu = Quu_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.condenseRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual_ref.lu -= dtau_ * s.beta;
  Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
  kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
                          + du_dq.transpose() * lu_condensed;
  kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
                          + du_dv.transpose() * lu_condensed; 
  kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
                          + du_da.transpose() * lu_condensed;
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isZero());
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Qfq().isZero());
  EXPECT_TRUE(kkt_matrix.Qfv().isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_residual.C().isZero());
  EXPECT_TRUE(kkt_matrix.Cq().isZero());
  EXPECT_TRUE(kkt_matrix.Cv().isZero());
  EXPECT_TRUE(kkt_matrix.Ca().isZero());
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
  SplitDirection d = SplitDirection::Random(robot);
  rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  std::cout << "du_dq" << std::endl;
  std::cout << du_dq << std::endl;
  std::cout << "du_dv" << std::endl;
  std::cout << du_dv << std::endl;
  std::cout << "du_da" << std::endl;
  std::cout << du_da << std::endl;
  std::cout << "Cq" << std::endl;
  std::cout << kkt_matrix.Cq() << std::endl;
  std::cout << "Cv" << std::endl;
  std::cout << kkt_matrix.Cv() << std::endl;
  std::cout << "Ca" << std::endl;
  std::cout << kkt_matrix.Ca() << std::endl;
  std::cout << "Cf" << std::endl;
  std::cout << kkt_matrix.Cf() << std::endl;
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
  const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimq(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimq(), robot.dimv());
  rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
  EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq));
  EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv));
}


TEST_F(RobotDynamicsTest, computeViolationL1NormFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {false};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.computeViolationL1Norm(robot, dtau_, s, kkt_residual), 
                   violation_ref);
}


TEST_F(RobotDynamicsTest, augmentRobotDynamicsFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {true};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(robot.dimf());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.augmentRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C());
  robot.computeBaumgarteDerivatives(dtau_, dtau_, kkt_matrix_ref.Cq(), 
                                    kkt_matrix_ref.Cv(), kkt_matrix_ref.Ca());
  kkt_residual_ref.lq() += dtau_ * du_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * du_dv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.la() += dtau_ * du_da.transpose() * s.beta 
                            + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dtau_ * du_df.transpose() * s.beta 
                            + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
  kkt_residual_ref.lu -= dtau_ * s.beta;
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C().topRows(robot.dimf())));
  const double violation_ref = kkt_residual_ref.C().lpNorm<1>() 
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
}


TEST_F(RobotDynamicsTest, condenseRobotDynamicsFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {true};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
  kkt_matrix.Quu = Quu_ref;
  kkt_matrix_ref.Quu = Quu_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(robot.dimf());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.condenseRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());  
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C());
  robot.computeBaumgarteDerivatives(dtau_, dtau_, kkt_matrix_ref.Cq(), 
                                    kkt_matrix_ref.Cv(), kkt_matrix_ref.Ca());
  kkt_residual_ref.lu -= dtau_ * s.beta;
  Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
  kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
                          + du_dq.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
                          + du_dv.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
                          + du_da.transpose() * lu_condensed 
                          + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
  kkt_residual_ref.lf() = dtau_ * du_df.transpose() * s.beta 
                          + du_df.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
  EXPECT_TRUE(kkt_matrix.Cf().isApprox(kkt_matrix_ref.Cf()));
  SplitDirection d = SplitDirection::Random(robot);
  rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  du_ref += du_df.leftCols(robot.dimf()) * d.df();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  std::cout << "du_dq" << std::endl;
  std::cout << du_dq << std::endl;
  std::cout << "du_dv" << std::endl;
  std::cout << du_dv << std::endl;
  std::cout << "du_da" << std::endl;
  std::cout << du_da << std::endl;
  std::cout << "du_df" << std::endl;
  std::cout << du_df << std::endl;
  std::cout << "Cq" << std::endl;
  std::cout << kkt_matrix.Cq() << std::endl;
  std::cout << "Cv" << std::endl;
  std::cout << kkt_matrix.Cv() << std::endl;
  std::cout << "Ca" << std::endl;
  std::cout << kkt_matrix.Ca() << std::endl;
  std::cout << "Cf" << std::endl;
  std::cout << kkt_matrix.Cf() << std::endl;
  const double violation_ref = kkt_residual_ref.C().lpNorm<1>() 
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
  const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
  EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq+du_df*df_dq));
  EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv+du_df*df_dv));
}


TEST_F(RobotDynamicsTest, computeViolationL1NormFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {true};
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C());
  kkt_residual_ref.u_res -= s.u;
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>()
                                + kkt_residual_ref.C().head(robot.dimf()).lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.computeViolationL1Norm(robot, dtau_, s, kkt_residual), 
                   violation_ref);
}


TEST_F(RobotDynamicsTest, augmentRobotDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(false);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  rd.augmentRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cu 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.dimf(), robot.dimv());  
  Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual_ref.lq() += dtau_ * du_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * du_dv.transpose() * s.beta; 
  kkt_residual_ref.la() += dtau_ * du_da.transpose() * s.beta;
  kkt_residual_ref.lu += - dtau_ * s.beta + Cu.transpose() * s.mu_stack();
  kkt_residual_ref.C().head(robot.dim_passive()) = dtau_ * s.u.head(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  const double violation_ref = dtau_ * s.u.head(6).lpNorm<1>()
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
}


TEST_F(RobotDynamicsTest, condenseRobotDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(false);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
  kkt_matrix.Quu = Quu_ref;
  kkt_matrix_ref.Quu = Quu_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.condenseRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());  
  Eigen::MatrixXd Cu 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.dimf(), robot.dimv());  
  Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual_ref.C().head(robot.dim_passive()) 
      = dtau_ * s.u.head(robot.dim_passive());
  kkt_residual_ref.lu -= dtau_ * s.beta;
  kkt_residual_ref.lu += Cu.transpose() * s.mu_stack();
  Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
  kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
                          + du_dq.transpose() * lu_condensed
                          + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
                          + du_dv.transpose() * lu_condensed
                          + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
                          + du_da.transpose() * lu_condensed
                          + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.C() += Cu * kkt_residual_ref.u_res;
  kkt_matrix_ref.Cq() += Cu * du_dq;
  kkt_matrix_ref.Cv() += Cu * du_dv;
  kkt_matrix_ref.Ca() += Cu * du_da;
  kkt_matrix_ref.Cf() += Cu * du_df;
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isZero());
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Qfq().isZero());
  EXPECT_TRUE(kkt_matrix.Qfv().isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
  SplitDirection d = SplitDirection::Random(robot);
  rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
  dbeta_ref += Cu.transpose() * d.dmu() / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  std::cout << "du_dq" << std::endl;
  std::cout << du_dq << std::endl;
  std::cout << "du_dv" << std::endl;
  std::cout << du_dv << std::endl;
  std::cout << "du_da" << std::endl;
  std::cout << du_da << std::endl;
  std::cout << "du_df" << std::endl;
  std::cout << du_df << std::endl;
  std::cout << "Cq" << std::endl;
  std::cout << kkt_matrix.Cq() << std::endl;
  std::cout << "Cv" << std::endl;
  std::cout << kkt_matrix.Cv() << std::endl;
  std::cout << "Ca" << std::endl;
  std::cout << kkt_matrix.Ca() << std::endl;
  std::cout << "Cf" << std::endl;
  std::cout << kkt_matrix.Cf() << std::endl;
  const double violation_ref = dtau_ * s.u.head(6).lpNorm<1>()
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
  const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
  EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq));
  EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv));
}


TEST_F(RobotDynamicsTest, computeViolationL1NormFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(false);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>()
                                + dtau_ * s.u.head(robot.dim_passive()).lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.computeViolationL1Norm(robot, dtau_, s, kkt_residual), 
                   violation_ref);
}


TEST_F(RobotDynamicsTest, augmentRobotDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(robot.dimf());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.augmentRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());  
  Eigen::MatrixXd Cu 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.dimf(), robot.dimv());  
  Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C().tail(robot.dimf()));
  robot.computeBaumgarteDerivatives(dtau_, dtau_, kkt_matrix_ref.Cq().bottomRows(robot.dimf()), 
                                    kkt_matrix_ref.Cv().bottomRows(robot.dimf()), 
                                    kkt_matrix_ref.Ca().bottomRows(robot.dimf()));
  kkt_residual_ref.lq() += dtau_ * du_dq.transpose() * s.beta 
                            + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * du_dv.transpose() * s.beta 
                            + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.la() += dtau_ * du_da.transpose() * s.beta 
                            + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dtau_ * du_df.transpose() * s.beta 
                            + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
  kkt_residual_ref.lu += - dtau_ * s.beta + Cu.transpose() * s.mu_stack();
  kkt_residual_ref.C().head(robot.dim_passive()) = dtau_ * s.u.head(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  const double violation_ref = kkt_residual_ref.C().tail(robot.dimf()).lpNorm<1>()
                                + dtau_ * s.u.head(6).lpNorm<1>()
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
}


TEST_F(RobotDynamicsTest, condenseRobotDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(robot);
  const Eigen::MatrixXd Quu_ref = Eigen::VectorXd::Random(robot.dimv()).asDiagonal();
  kkt_matrix.Quu = Quu_ref;
  kkt_matrix_ref.Quu = Quu_ref;
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(robot.dimf());
  kkt_residual.lu = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la() = kkt_residual.la();
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu = kkt_residual.lu;
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  rd.condenseRobotDynamics(robot, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());  
  Eigen::MatrixXd Cu 
      = Eigen::MatrixXd::Zero(robot.dim_passive()+robot.dimf(), robot.dimv());  
  Cu.topLeftCorner(robot.dim_passive(), robot.dim_passive()).diagonal().fill(dtau_);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  kkt_residual_ref.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C().tail(robot.dimf()));
  kkt_residual_ref.C().head(robot.dim_passive()) 
      = dtau_ * s.u.head(robot.dim_passive());
  robot.computeBaumgarteDerivatives(dtau_, dtau_, kkt_matrix_ref.Cq().bottomRows(robot.dimf()), 
                                    kkt_matrix_ref.Cv().bottomRows(robot.dimf()), 
                                    kkt_matrix_ref.Ca().bottomRows(robot.dimf()));
  kkt_residual_ref.lu -= dtau_ * s.beta;
  kkt_residual_ref.lu += Cu.transpose() * s.mu_stack();
  Eigen::VectorXd lu_condensed = kkt_residual_ref.lu + kkt_matrix_ref.Quu * kkt_residual_ref.u_res; 
  kkt_residual_ref.lq() = dtau_ * du_dq.transpose() * s.beta 
                          + du_dq.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() = dtau_ * du_dv.transpose() * s.beta 
                          + du_dv.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cv().transpose() * s.mu_stack();
  kkt_residual_ref.la() = dtau_ * du_da.transpose() * s.beta 
                          + du_da.transpose() * lu_condensed 
                          + kkt_matrix_ref.Ca().transpose() * s.mu_stack();
  kkt_residual_ref.lf() = dtau_ * du_df.transpose() * s.beta 
                          + du_df.transpose() * lu_condensed 
                          + kkt_matrix_ref.Cf().transpose() * s.mu_stack();
  kkt_residual_ref.C() += Cu * kkt_residual_ref.u_res;
  kkt_matrix_ref.Cq() += Cu * du_dq;
  kkt_matrix_ref.Cv() += Cu * du_dv;
  kkt_matrix_ref.Ca() += Cu * du_da;
  kkt_matrix_ref.Cf() += Cu * du_df;
  EXPECT_TRUE(kkt_residual.u_res.isApprox(kkt_residual_ref.u_res));
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la().isApprox(kkt_residual_ref.la()));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu.isApprox(kkt_residual_ref.lu));
  EXPECT_TRUE(kkt_matrix.Qaa().isApprox(du_da.transpose()*Quu_ref*du_da));
  EXPECT_TRUE(kkt_matrix.Qaf().isApprox(du_da.transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qaq().isApprox(du_da.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qav().isApprox(du_da.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qff().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_df.leftCols(robot.dimf())));
  EXPECT_TRUE(kkt_matrix.Qfq().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qfv().isApprox(du_df.leftCols(robot.dimf()).transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(du_dq.transpose()*Quu_ref*du_dq));
  EXPECT_TRUE(kkt_matrix.Qqv().isApprox(du_dq.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_matrix.Qvv().isApprox(du_dv.transpose()*Quu_ref*du_dv));
  EXPECT_TRUE(kkt_residual.C().isApprox(kkt_residual_ref.C()));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(kkt_matrix_ref.Cq()));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(kkt_matrix_ref.Cv()));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(kkt_matrix_ref.Ca()));
  EXPECT_TRUE(kkt_matrix.Cf().isApprox(kkt_matrix_ref.Cf()));
  SplitDirection d = SplitDirection::Random(robot);
  rd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd du_ref = kkt_residual_ref.u_res;
  du_ref += du_dq * d.dq();
  du_ref += du_dv * d.dv();
  du_ref += du_da * d.da();
  if (robot.has_active_contacts()) {
    du_ref += du_df.leftCols(robot.dimf()) * d.df();
  }
  EXPECT_TRUE(du_ref.isApprox(d.du));
  Eigen::VectorXd dbeta_ref = (kkt_residual_ref.lu + Quu_ref * du_ref) / dtau_;
  dbeta_ref += Cu.transpose() * d.dmu() / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta));
  std::cout << "du_dq" << std::endl;
  std::cout << du_dq << std::endl;
  std::cout << "du_dv" << std::endl;
  std::cout << du_dv << std::endl;
  std::cout << "du_da" << std::endl;
  std::cout << du_da << std::endl;
  std::cout << "du_df" << std::endl;
  std::cout << du_df << std::endl;
  std::cout << "Cq" << std::endl;
  std::cout << kkt_matrix.Cq() << std::endl;
  std::cout << "Cv" << std::endl;
  std::cout << kkt_matrix.Cv() << std::endl;
  std::cout << "Ca" << std::endl;
  std::cout << kkt_matrix.Ca() << std::endl;
  std::cout << "Cf" << std::endl;
  std::cout << kkt_matrix.Cf() << std::endl;
  std::cout << "Cu" << std::endl;
  std::cout << Cu << std::endl;
  std::cout << "Cu.transpose()" << std::endl;
  std::cout << Cu.transpose() << std::endl;
  const double violation_ref = kkt_residual_ref.C().tail(robot.dimf()).lpNorm<1>()
                                + dtau_ * s.u.head(6).lpNorm<1>()
                                + dtau_ * kkt_residual_ref.u_res.lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.violationL1Norm(robot, dtau_, s, kkt_residual), violation_ref);
  const Eigen::MatrixXd da_dq = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd da_dv = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd df_dq = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  const Eigen::MatrixXd df_dv = Eigen::MatrixXd::Random(robot.dimf(), robot.dimv());
  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  rd.getStateFeedbackGain(da_dq, da_dv, df_dq, df_dv, Kuq, Kuv);
  EXPECT_TRUE(Kuq.isApprox(du_dq+du_da*da_dq+du_df*df_dq));
  EXPECT_TRUE(Kuv.isApprox(du_dv+du_da*da_dv+du_df*df_dv));
}


TEST_F(RobotDynamicsTest, computeViolationL1NormFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s = SplitSolution::Random(robot);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(robot);
  RobotDynamics rd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual_ref.u_res);
  robot.computeBaumgarteResidual(dtau_, dtau_, kkt_residual_ref.C().tail(robot.dimf()));
  kkt_residual_ref.u_res -= s.u;
  const double violation_ref = dtau_ * kkt_residual_ref.u_res.lpNorm<1>()
                                + dtau_ * s.u.head(robot.dim_passive()).lpNorm<1>()
                                + kkt_residual_ref.C().tail(robot.dimf()).lpNorm<1>();
  EXPECT_DOUBLE_EQ(rd.computeViolationL1Norm(robot, dtau_, s, kkt_residual), 
                   violation_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}