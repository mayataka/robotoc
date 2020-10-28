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
#include "idocp/ocp/contact_dynamics_data.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

class ContactDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


TEST_F(ContactDynamicsTest, linearizeInverseDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID_full().noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


TEST_F(ContactDynamicsTest, linearizeContactConstraintFixedBaseWitoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  EXPECT_TRUE(data.IDC().isZero());
  EXPECT_TRUE(data.dCda().isZero());
  EXPECT_TRUE(data.dIDCdqv().isZero());
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  kkt_residual_ref.lq() += dtau * data.dIDdq().transpose() * s.beta;
  kkt_residual_ref.lv() += dtau * data.dIDdv().transpose() * s.beta; 
  kkt_residual_ref.la += dtau * data.dIDda.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau * s.beta;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm_ref = dtau * data.ID().lpNorm<1>();
  const double squarednorm_ref = dtau * dtau * data.ID().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeInverseDynamicsFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID_full().noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


TEST_F(ContactDynamicsTest, linearizeContactConstraintFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  data.setContactStatus(contact_status);
  data_ref.setContactStatus(contact_status);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau, data_ref.C());
  robot.computeBaumgarteDerivatives(contact_status, dtau, data_ref.dCdq(), 
                                    data_ref.dCdv(), data_ref.dCda());
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
  EXPECT_TRUE(data.dCda().isApprox(data_ref.dCda()));
  EXPECT_TRUE(data.dIDCdqv().isApprox(data_ref.dIDCdqv()));
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  Eigen::MatrixXd dIDdf = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dIDdf);
  kkt_residual_ref.lq() += dtau * data.dIDdq().transpose() * s.beta + dtau * data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau * data.dIDdv().transpose() * s.beta + dtau * data.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau * data.dIDda.transpose() * s.beta + dtau * data.dCda().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dtau * dIDdf.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau * s.beta;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm_ref = dtau * data.IDC().lpNorm<1>();
  const double squarednorm_ref = dtau * dtau * data.IDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeInverseDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false, false, false, false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID_full().head(robot.dim_passive()).noalias() -= s.u_passive;
  data_ref.ID_full().tail(robot.dimu()).noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


TEST_F(ContactDynamicsTest, linearizeContactConstraintFloatingBaseWitoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false, false, false, false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  EXPECT_TRUE(data.IDC().isZero());
  EXPECT_TRUE(data.dCda().isZero());
  EXPECT_TRUE(data.dIDCdqv().isZero());
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false, false, false, false};
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(robot.dimv());
  lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  ContactDynamics cd(robot);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  kkt_residual_ref.lq() += dtau * data.dIDdq().transpose() * s.beta;
  kkt_residual_ref.lv() += dtau * data.dIDdv().transpose() * s.beta; 
  kkt_residual_ref.la += dtau * data.dIDda.transpose() * s.beta;
  lu_full_ref -= dtau * s.beta;
  lu_full_ref.head(robot.dim_passive()) += dtau * s.nu_passive;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(lu_full_ref.head(robot.dim_passive())));
  EXPECT_TRUE(kkt_residual.lu().isApprox(lu_full_ref.tail(robot.dimu())));
  const double l1norm_ref = dtau * (data.ID_full().lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau * dtau * (data.ID_full().squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeInverseDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID_full().head(robot.dim_passive()).noalias() -= s.u_passive;
  data_ref.ID_full().tail(robot.dimu()).noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


TEST_F(ContactDynamicsTest, linearizeContactConstraintFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  robot.updateKinematics(s.q, s.v, s.a);
  data_ref.setContactStatus(contact_status);
  robot.computeBaumgarteResidual(contact_status, dtau, data_ref.C());
  robot.computeBaumgarteDerivatives(contact_status, dtau, data_ref.dCdq(), 
                                    data_ref.dCdv(), data_ref.dCda());
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
  EXPECT_TRUE(data.dCda().isApprox(data_ref.dCda()));
  EXPECT_TRUE(data.dIDCdqv().isApprox(data_ref.dIDCdqv()));
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(robot.dimv());
  lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  Eigen::MatrixXd dIDdf = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dIDdf);
  kkt_residual_ref.lq() += dtau * data.dIDdq().transpose() * s.beta + dtau * data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau * data.dIDdv().transpose() * s.beta + dtau * data.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau * data.dIDda.transpose() * s.beta + dtau * data.dCda().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dtau * dIDdf.transpose() * s.beta;
  lu_full_ref -= dtau * s.beta;
  lu_full_ref.head(robot.dim_passive()) += dtau * s.nu_passive;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(lu_full_ref.head(robot.dim_passive())));
  EXPECT_TRUE(kkt_residual.lu().isApprox(lu_full_ref.tail(robot.dimu())));
  const double l1norm_ref = dtau * (data.IDC().lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau * dtau * (data.IDC().squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condensingFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  ContactDynamicsData data_ref = data;
  ContactDynamics::condensing(robot, dtau, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC() = data_ref.MJtJinv() * data_ref.IDC();
  data_ref.Qafqv() = - kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  data_ref.Qafu() = kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * IO_mat;
  data_ref.laf().head(robot.dimv()) = kkt_residual_ref.la;
  data_ref.laf().tail(contact_status.dimf()) = kkt_residual_ref.lf();
  data_ref.laf() -= kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxu() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu();
  kkt_matrix_ref.Quu() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu();
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  kkt_residual_ref.lu() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvq() = - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) 
                          - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = (OOIO_mat * data_ref.MJtJinv() * IO_mat).bottomRows(robot.dimv());
  kkt_residual_ref.Fx() -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());
  EXPECT_TRUE(data_ref.MJtJinv_dIDCdqv().isApprox(data.MJtJinv_dIDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_IDC().isApprox(data.MJtJinv_IDC()));
  EXPECT_TRUE(data_ref.Qafqv().isApprox(data.Qafqv()));
  EXPECT_TRUE(data_ref.Qafu().isApprox(data.Qafu()));
  EXPECT_TRUE(data_ref.laf().isApprox(data.laf()));
  EXPECT_TRUE(kkt_matrix_ref.Qxx().isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix_ref.Qxu().isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(kkt_matrix_ref.Quu().isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_residual_ref.lx().isApprox(kkt_residual.lx()));
  EXPECT_TRUE(kkt_residual_ref.lu().isApprox(kkt_residual.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isZero());
  EXPECT_TRUE(kkt_matrix_ref.Fqq().isApprox(kkt_matrix.Fqq()));
  EXPECT_TRUE(kkt_matrix_ref.Fqv().isApprox(kkt_matrix.Fqv()));
  EXPECT_TRUE(kkt_matrix_ref.Fqu().isApprox(kkt_matrix.Fqu()));
  EXPECT_TRUE(kkt_matrix_ref.Fvq().isApprox(kkt_matrix.Fvq()));
  EXPECT_TRUE(kkt_matrix_ref.Fvv().isApprox(kkt_matrix.Fvv()));
  EXPECT_TRUE(kkt_matrix_ref.Fvu().isApprox(kkt_matrix.Fvu()));
  EXPECT_TRUE(kkt_residual_ref.Fx().isApprox(kkt_residual.Fx()));
  std::cout << "IO_mat" << std::endl;
  std::cout << IO_mat << std::endl;
  std::cout << "OOIO_mat" << std::endl;
  std::cout << OOIO_mat << std::endl;
}


TEST_F(ContactDynamicsTest, condensingFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  ContactDynamicsData data_ref = data;
  ContactDynamics::condensing(robot, dtau, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC() = data_ref.MJtJinv() * data_ref.IDC();
  data_ref.Qafqv() = - kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  data_ref.Qafu() = kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * IO_mat;
  data_ref.laf().head(robot.dimv()) = kkt_residual_ref.la;
  data_ref.laf().tail(contact_status.dimf()) = kkt_residual_ref.lf();
  data_ref.laf() -= kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxu() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu();
  kkt_matrix_ref.Quu() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu();
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  kkt_residual_ref.lu() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvq() = - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) 
                          - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = (OOIO_mat * data_ref.MJtJinv() * IO_mat).bottomRows(robot.dimv());
  kkt_residual_ref.Fx() -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());
  EXPECT_TRUE(data_ref.MJtJinv_dIDCdqv().isApprox(data.MJtJinv_dIDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_IDC().isApprox(data.MJtJinv_IDC()));
  EXPECT_TRUE(data_ref.Qafqv().isApprox(data.Qafqv()));
  EXPECT_TRUE(data_ref.Qafu().isApprox(data.Qafu()));
  EXPECT_TRUE(data_ref.laf().isApprox(data.laf()));
  EXPECT_TRUE(kkt_matrix_ref.Qxx().isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix_ref.Qxu().isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(kkt_matrix_ref.Quu().isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_residual_ref.lx().isApprox(kkt_residual.lx()));
  EXPECT_TRUE(kkt_residual_ref.lu().isApprox(kkt_residual.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isZero());
  EXPECT_TRUE(kkt_matrix_ref.Fqq().isApprox(kkt_matrix.Fqq()));
  EXPECT_TRUE(kkt_matrix_ref.Fqv().isApprox(kkt_matrix.Fqv()));
  EXPECT_TRUE(kkt_matrix_ref.Fqu().isApprox(kkt_matrix.Fqu()));
  EXPECT_TRUE(kkt_matrix_ref.Fvq().isApprox(kkt_matrix.Fvq()));
  EXPECT_TRUE(kkt_matrix_ref.Fvv().isApprox(kkt_matrix.Fvv()));
  EXPECT_TRUE(kkt_matrix_ref.Fvu().isApprox(kkt_matrix.Fvu()));
  EXPECT_TRUE(kkt_residual_ref.Fx().isApprox(kkt_residual.Fx()));
  std::cout << "IO_mat" << std::endl;
  std::cout << IO_mat << std::endl;
  std::cout << "OOIO_mat" << std::endl;
  std::cout << OOIO_mat << std::endl;
}


TEST_F(ContactDynamicsTest, condensingFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false, false, false, false};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  ContactDynamicsData data_ref = data;
  ContactDynamics::condensing(robot, dtau, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC() = data_ref.MJtJinv() * data_ref.IDC();
  data_ref.Qafqv() = - kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  data_ref.Qafu_full() = kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * IO_mat;
  data_ref.laf().head(robot.dimv()) = kkt_residual_ref.la;
  data_ref.laf().tail(contact_status.dimf()) = kkt_residual_ref.lf();
  data_ref.laf() -= kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxu_full() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu_full();
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu_full();
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(robot.dimv());
  lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  lu_full_ref += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  kkt_residual_ref.lu_passive = lu_full_ref.head(robot.dim_passive());
  kkt_residual_ref.lu() = lu_full_ref.tail(robot.dimu());
  kkt_residual_ref.lu() -= kkt_matrix_ref.Quu_passive_bottomLeft() * data_ref.u_passive;
  kkt_residual_ref.lx() -= kkt_matrix_ref.Qxu_passive() * data_ref.u_passive;
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvq() = - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) 
                          - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomRightCorner(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data_ref.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fvu() = (Fxu_full.rightCols(robot.dimu())).bottomRows(robot.dimv());
  kkt_residual_ref.Fx() -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());
  kkt_residual_ref.Fx() -= Fxu_full.leftCols(robot.dim_passive()) * data_ref.u_passive;
  EXPECT_TRUE(data_ref.MJtJinv_dIDCdqv().isApprox(data.MJtJinv_dIDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_IDC().isApprox(data.MJtJinv_IDC()));
  EXPECT_TRUE(data_ref.Qafqv().isApprox(data.Qafqv()));
  EXPECT_TRUE(data_ref.Qafu().isApprox(data.Qafu()));
  EXPECT_TRUE(data_ref.laf().isApprox(data.laf()));
  EXPECT_TRUE(kkt_matrix_ref.Qxx().isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix_ref.Qxu().isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(kkt_matrix_ref.Quu().isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_residual_ref.lx().isApprox(kkt_residual.lx()));
  EXPECT_TRUE(kkt_residual_ref.lu().isApprox(kkt_residual.lu()));
  EXPECT_TRUE(kkt_residual_ref.lu_passive.isApprox(kkt_residual.lu_passive));
  EXPECT_TRUE(kkt_matrix_ref.Fqq().isApprox(kkt_matrix.Fqq()));
  EXPECT_TRUE(kkt_matrix_ref.Fqv().isApprox(kkt_matrix.Fqv()));
  EXPECT_TRUE(kkt_matrix_ref.Fqu().isApprox(kkt_matrix.Fqu()));
  EXPECT_TRUE(kkt_matrix_ref.Fvq().isApprox(kkt_matrix.Fvq()));
  EXPECT_TRUE(kkt_matrix_ref.Fvv().isApprox(kkt_matrix.Fvv()));
  EXPECT_TRUE(kkt_matrix_ref.Fvu().isApprox(kkt_matrix.Fvu()));
  EXPECT_TRUE(kkt_residual_ref.Fx().isApprox(kkt_residual.Fx()));
  std::cout << "IO_mat" << std::endl;
  std::cout << IO_mat << std::endl;
  std::cout << "OOIO_mat" << std::endl;
  std::cout << OOIO_mat << std::endl;
}



TEST_F(ContactDynamicsTest, condensingFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  ContactDynamicsData data_ref = data;
  ContactDynamics::condensing(robot, dtau, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC() = data_ref.MJtJinv() * data_ref.IDC();
  data_ref.Qafqv() = - kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  data_ref.Qafu_full() = kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * IO_mat;
  data_ref.laf().head(robot.dimv()) = kkt_residual_ref.la;
  data_ref.laf().tail(contact_status.dimf()) = kkt_residual_ref.lf();
  data_ref.laf() -= kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxu_full() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu_full();
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu_full();
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(robot.dimv());
  lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  lu_full_ref += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  kkt_residual_ref.lu_passive = lu_full_ref.head(robot.dim_passive());
  kkt_residual_ref.lu() = lu_full_ref.tail(robot.dimu());
  kkt_residual_ref.lu() -= kkt_matrix_ref.Quu_passive_bottomLeft() * data_ref.u_passive;
  kkt_residual_ref.lx() -= kkt_matrix_ref.Qxu_passive() * data_ref.u_passive;
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvq() = - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) 
                          - (OOIO_mat * data_ref.MJtJinv_dIDCdqv()).bottomRightCorner(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data_ref.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fvu() = (Fxu_full.rightCols(robot.dimu())).bottomRows(robot.dimv());
  kkt_residual_ref.Fx() -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());
  kkt_residual_ref.Fx() -= Fxu_full.leftCols(robot.dim_passive()) * data_ref.u_passive;
  EXPECT_TRUE(data_ref.MJtJinv_dIDCdqv().isApprox(data.MJtJinv_dIDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_IDC().isApprox(data.MJtJinv_IDC()));
  EXPECT_TRUE(data_ref.Qafqv().isApprox(data.Qafqv()));
  EXPECT_TRUE(data_ref.Qafu().isApprox(data.Qafu()));
  EXPECT_TRUE(data_ref.laf().isApprox(data.laf()));
  EXPECT_TRUE(kkt_matrix_ref.Qxx().isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix_ref.Qxu().isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(kkt_matrix_ref.Quu().isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_residual_ref.lx().isApprox(kkt_residual.lx()));
  EXPECT_TRUE(kkt_residual_ref.lu().isApprox(kkt_residual.lu()));
  EXPECT_TRUE(kkt_residual_ref.lu_passive.isApprox(kkt_residual.lu_passive));
  EXPECT_TRUE(kkt_matrix_ref.Fqq().isApprox(kkt_matrix.Fqq()));
  EXPECT_TRUE(kkt_matrix_ref.Fqv().isApprox(kkt_matrix.Fqv()));
  EXPECT_TRUE(kkt_matrix_ref.Fqu().isApprox(kkt_matrix.Fqu()));
  EXPECT_TRUE(kkt_matrix_ref.Fvq().isApprox(kkt_matrix.Fvq()));
  EXPECT_TRUE(kkt_matrix_ref.Fvv().isApprox(kkt_matrix.Fvv()));
  EXPECT_TRUE(kkt_matrix_ref.Fvu().isApprox(kkt_matrix.Fvu()));
  EXPECT_TRUE(kkt_residual_ref.Fx().isApprox(kkt_residual.Fx()));
  std::cout << "IO_mat" << std::endl;
  std::cout << IO_mat << std::endl;
  std::cout << "OOIO_mat" << std::endl;
  std::cout << OOIO_mat << std::endl;
}


TEST_F(ContactDynamicsTest, expansionFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  const ContactDynamicsData data_ref = data;
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*robot.dimv());
  ContactDynamics::expansion(robot, dtau, data, kkt_matrix, kkt_residual, dlmdgmm_next.tail(robot.dimv()), d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IO_mat.topRows(robot.dimv()).setIdentity();
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  d_ref.daf() = - data_ref.MJtJinv_dIDCdqv() * d_ref.dx() + data_ref.MJtJinv() * IO_mat * d_ref.du()
                  - data_ref.MJtJinv_IDC();
  d_ref.df().array() *= -1;
  d_ref.dbetamu() = - data_ref.MJtJinv() * data_ref.Qafqv() * d_ref.dx()
                    - data_ref.MJtJinv() * data_ref.Qafu() * d_ref.du()
                    - data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next
                    - data_ref.MJtJinv() * data_ref.laf();
  d_ref.dbetamu() /= dtau;
  EXPECT_TRUE(d_ref.da().isApprox(d.da()));
  EXPECT_TRUE(d_ref.df().isApprox(d.df()));
  EXPECT_TRUE(d_ref.dbeta().isApprox(d.dbeta()));
  EXPECT_TRUE(d_ref.dmu().isApprox(d.dmu()));
}


TEST_F(ContactDynamicsTest, expansionFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  const ContactDynamicsData data_ref = data;
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*robot.dimv());
  ContactDynamics::expansion(robot, dtau, data, kkt_matrix, kkt_residual, dlmdgmm_next.tail(robot.dimv()), d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IO_mat.topRows(robot.dimv()).setIdentity();
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  d_ref.daf() = - data_ref.MJtJinv_dIDCdqv() * d_ref.dx() + data_ref.MJtJinv() * IO_mat * d_ref.du()
                  - data_ref.MJtJinv_IDC();
  d_ref.df().array() *= -1;
  d_ref.dbetamu() = - data_ref.MJtJinv() * data_ref.Qafqv() * d_ref.dx()
                    - data_ref.MJtJinv() * data_ref.Qafu() * d_ref.du()
                    - data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next
                    - data_ref.MJtJinv() * data_ref.laf();
  d_ref.dbetamu() /= dtau;
  EXPECT_TRUE(d_ref.da().isApprox(d.da()));
  EXPECT_TRUE(d_ref.df().isApprox(d.df()));
  EXPECT_TRUE(d_ref.dbeta().isApprox(d.dbeta()));
  EXPECT_TRUE(d_ref.dmu().isApprox(d.dmu()));
}


TEST_F(ContactDynamicsTest, expansionFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false, false, false, false};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.lu_passive = Eigen::VectorXd::Random(robot.dim_passive());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  const ContactDynamicsData data_ref = data;
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*robot.dimv());
  ContactDynamics::expansion(robot, dtau, data, kkt_matrix, kkt_residual, dlmdgmm_next.tail(robot.dimv()), d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IO_mat.topRows(robot.dimv()).setIdentity();
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  d_ref.du_passive = - data_ref.u_passive;
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(robot.dimv());
  du_full.head(robot.dim_passive()) = d_ref.du_passive;
  du_full.tail(robot.dimu()) = d_ref.du();
  d_ref.dnu_passive = - kkt_residual_ref.lu_passive 
                      - kkt_matrix_ref.Quu_full().topRows(robot.dim_passive()) * du_full
                      - (kkt_matrix_ref.Qxu_full().transpose()).topRows(robot.dim_passive()) * d_ref.dx()
                      - (IO_mat.transpose()).topRows(robot.dim_passive()) * data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next;
  d_ref.dnu_passive /= dtau;
  d_ref.daf() = - data_ref.MJtJinv_dIDCdqv() * d_ref.dx() + data_ref.MJtJinv() * IO_mat * du_full
                  - data_ref.MJtJinv_IDC();
  d_ref.df().array() *= -1;
  d_ref.dbetamu() = - data_ref.MJtJinv() * data_ref.Qafqv() * d_ref.dx()
                    - data_ref.MJtJinv() * data_ref.Qafu_full() * du_full
                    - data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next
                    - data_ref.MJtJinv() * data_ref.laf();
  d_ref.dbetamu() /= dtau;
  EXPECT_TRUE(d_ref.du_passive.isApprox(d.du_passive));
  EXPECT_TRUE(d_ref.dnu_passive.isApprox(d.dnu_passive));
  EXPECT_TRUE(d_ref.da().isApprox(d.da()));
  EXPECT_TRUE(d_ref.df().isApprox(d.df()));
  EXPECT_TRUE(d_ref.dbeta().isApprox(d.dbeta()));
  EXPECT_TRUE(d_ref.dmu().isApprox(d.dmu()));
}


TEST_F(ContactDynamicsTest, expansionFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual.lu_passive = Eigen::VectorXd::Random(robot.dim_passive());
  kkt_residual.Fx() = Eigen::VectorXd::Random(2*robot.dimv());
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  KKTResidual kkt_residual_ref = kkt_residual;
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  data.u_passive.setRandom();
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  const ContactDynamicsData data_ref = data;
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*robot.dimv());
  ContactDynamics::expansion(robot, dtau, data, kkt_matrix, kkt_residual, dlmdgmm_next.tail(robot.dimv()), d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IO_mat.topRows(robot.dimv()).setIdentity();
  OOIO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  d_ref.du_passive = - data_ref.u_passive;
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(robot.dimv());
  du_full.head(robot.dim_passive()) = d_ref.du_passive;
  du_full.tail(robot.dimu()) = d_ref.du();
  d_ref.dnu_passive = - kkt_residual_ref.lu_passive 
                      - kkt_matrix_ref.Quu_full().topRows(robot.dim_passive()) * du_full
                      - (kkt_matrix_ref.Qxu_full().transpose()).topRows(robot.dim_passive()) * d_ref.dx()
                      - (IO_mat.transpose()).topRows(robot.dim_passive()) * data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next;
  d_ref.dnu_passive /= dtau;
  d_ref.daf() = - data_ref.MJtJinv_dIDCdqv() * d_ref.dx() + data_ref.MJtJinv() * IO_mat * du_full
                  - data_ref.MJtJinv_IDC();
  d_ref.df().array() *= -1;
  d_ref.dbetamu() = - data_ref.MJtJinv() * data_ref.Qafqv() * d_ref.dx()
                    - data_ref.MJtJinv() * data_ref.Qafu_full() * du_full
                    - data_ref.MJtJinv() * OOIO_mat.transpose() * dlmdgmm_next
                    - data_ref.MJtJinv() * data_ref.laf();
  d_ref.dbetamu() /= dtau;
  EXPECT_TRUE(d_ref.du_passive.isApprox(d.du_passive));
  EXPECT_TRUE(d_ref.dnu_passive.isApprox(d.dnu_passive));
  EXPECT_TRUE(d_ref.da().isApprox(d.da()));
  EXPECT_TRUE(d_ref.df().isApprox(d.df()));
  EXPECT_TRUE(d_ref.dbeta().isApprox(d.dbeta()));
  EXPECT_TRUE(d_ref.dmu().isApprox(d.dmu()));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}