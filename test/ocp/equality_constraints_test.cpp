#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/equality_constraints.hpp"


namespace idocp {

class EqualityConstraintsTest : public ::testing::Test {
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


TEST_F(EqualityConstraintsTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(fixed_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  // std::vector<bool> contact_status = {rnd()%2==0};
  std::vector<bool> contact_status = {true};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau_, s, 
                                                    kkt_matrix, kkt_residual);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd la = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.dimf()+robot.dim_passive());
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimf());
  robot.computeBaumgarteResidual(dtau_, C);
  robot.computeBaumgarteDerivatives(dtau_, Cq, Cv, Ca);
  lq += Cq.transpose() * s.mu_active();
  lv += Cv.transpose() * s.mu_active();
  la += Ca.transpose() * s.mu_active();
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv));
  EXPECT_TRUE(kkt_residual.la().isApprox(la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(Cq));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(Cv));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(Ca));
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
}



TEST_F(EqualityConstraintsTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
  Robot robot(floating_base_urdf_, contact_frames, baum_a, baum_b);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  robot.generateFeasibleConfiguration(s.q);
  s.v = Eigen::VectorXd::Random(robot.dimv());
  s.a = Eigen::VectorXd::Random(robot.dimv());
  s.f = Eigen::VectorXd::Random(robot.max_dimf());
  s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
  s.lmd = Eigen::VectorXd::Random(robot.dimv());
  s.gmm = Eigen::VectorXd::Random(robot.dimv());
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(robot);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau_, s, 
                                                    kkt_matrix, kkt_residual);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd la = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(robot.dimf()+robot.dim_passive());
  Eigen::MatrixXd Cq = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cv = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Ca = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimv());
  Eigen::MatrixXd Cf = Eigen::MatrixXd::Zero(robot.dimf()+robot.dim_passive(), robot.dimf());
  robot.computeBaumgarteResidual(dtau_, C);
  robot.computeBaumgarteDerivatives(dtau_, Cq, Cv, Ca);
  C.tail(robot.dim_passive()) += s.u.head(robot.dim_passive());
  lq += Cq.transpose() * s.mu_active();
  lv += Cv.transpose() * s.mu_active();
  la += Ca.transpose() * s.mu_active();
  lu.head(robot.dim_passive()) += dtau_ * s.mu_active().tail(robot.dim_passive());
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv));
  EXPECT_TRUE(kkt_residual.la().isApprox(la));
  EXPECT_TRUE(kkt_residual.lu.isApprox(lu));
  EXPECT_TRUE(kkt_matrix.Cq().isApprox(Cq));
  EXPECT_TRUE(kkt_matrix.Cv().isApprox(Cv));
  EXPECT_TRUE(kkt_matrix.Ca().isApprox(Ca));
  EXPECT_TRUE(kkt_matrix.Cf().isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}