#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"


namespace idocp {

class SplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(SplitSolutionTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu = Eigen::VectorXd::Random(robot.max_dimf()+robot.dim_passive());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd f = Eigen::VectorXd::Random(robot.max_dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu = mu;
  s.a = a;
  s.f = f;
  s.q = q;
  s.v = v;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu.isApprox(mu));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.f.isApprox(f));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.mu_active().isApprox(mu.head(robot.dimf()+robot.dim_passive())));
  EXPECT_TRUE(s.f_active().isApprox(f.head(robot.dimf())));
  s.mu_active().setZero();
  s.f_active().setZero();
  EXPECT_TRUE(s.mu_active().isZero());
  EXPECT_TRUE(s.f_active().isZero());
  if (robot.dimf() < robot.max_dimf()) {
    EXPECT_FALSE(s.mu.isZero());
    EXPECT_FALSE(s.f.isZero());
  }
  s.mu.setZero();
  s.f.setZero();
  const Eigen::VectorXd mu_active = Eigen::VectorXd::Random(robot.dimf()+robot.dim_passive());
  const Eigen::VectorXd f_active = Eigen::VectorXd::Random(robot.dimf());
  s.mu_active() = mu_active;
  s.f_active() = f_active;
  EXPECT_TRUE(s.mu_active().isApprox(mu_active));
  EXPECT_TRUE(s.f_active().isApprox(f_active));
}


TEST_F(SplitSolutionTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu = Eigen::VectorXd::Random(robot.max_dimf()+robot.dim_passive());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd f = Eigen::VectorXd::Random(robot.max_dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu = mu;
  s.a = a;
  s.f = f;
  s.q = q;
  s.v = v;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu.isApprox(mu));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.f.isApprox(f));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.mu_active().isApprox(mu.head(robot.dimf()+robot.dim_passive())));
  EXPECT_TRUE(s.f_active().isApprox(f.head(robot.dimf())));
  s.mu_active().setZero();
  s.f_active().setZero();
  EXPECT_TRUE(s.mu_active().isZero());
  EXPECT_TRUE(s.f_active().isZero());
  if (robot.dimf() < robot.max_dimf()) {
    EXPECT_FALSE(s.mu.isZero());
    EXPECT_FALSE(s.f.isZero());
  }
  s.mu.setZero();
  s.f.setZero();
  const Eigen::VectorXd mu_active = Eigen::VectorXd::Random(robot.dimf()+robot.dim_passive());
  const Eigen::VectorXd f_active = Eigen::VectorXd::Random(robot.dimf());
  s.mu_active() = mu_active;
  s.f_active() = f_active;
  EXPECT_TRUE(s.mu_active().isApprox(mu_active));
  EXPECT_TRUE(s.f_active().isApprox(f_active));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}