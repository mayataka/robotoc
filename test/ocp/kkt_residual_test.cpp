#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/kkt_residual.hpp"


namespace idocp {

class KKTResidualTest : public ::testing::Test {
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


TEST_F(KKTResidualTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  KKTComposition composition(robot);
  composition.set(robot);
  KKTResidual residual(robot);
  residual.setContactStatus(robot);
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(composition.Fq_size());
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(composition.Fv_size());
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(composition.C_size());
  const Eigen::VectorXd Qa_res = Eigen::VectorXd::Random(composition.Qa_size());
  const Eigen::VectorXd Qf_res = Eigen::VectorXd::Random(composition.Qf_size());
  const Eigen::VectorXd Qq_res = Eigen::VectorXd::Random(composition.Qq_size());
  const Eigen::VectorXd Qv_res = Eigen::VectorXd::Random(composition.Qv_size());
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.la() = Qa_res;
  residual.lf() = Qf_res;
  residual.lq() = Qq_res;
  residual.lv() = Qv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Fq_begin(), composition.Fq_size()).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Fv_begin(), composition.Fv_size()).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.C_begin(), composition.C_size()).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qa_begin(), composition.Qa_size()).isApprox(Qa_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qf_begin(), composition.Qf_size()).isApprox(Qf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qq_begin(), composition.Qq_size()).isApprox(Qq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qv_begin(), composition.Qv_size()).isApprox(Qv_res));
  EXPECT_TRUE(residual.lx().head(composition.Qq_size()).isApprox(Qq_res));
  EXPECT_TRUE(residual.lx().tail(composition.Qv_size()).isApprox(Qv_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimKKT(), composition.dimKKT());
  EXPECT_EQ(residual.max_dimKKT(), composition.max_dimKKT());
}


TEST_F(KKTResidualTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  KKTComposition composition(robot);
  composition.set(robot);
  KKTResidual residual(robot);
  residual.setContactStatus(robot);
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(composition.Fq_size());
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(composition.Fv_size());
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(composition.C_size());
  const Eigen::VectorXd Qa_res = Eigen::VectorXd::Random(composition.Qa_size());
  const Eigen::VectorXd Qf_res = Eigen::VectorXd::Random(composition.Qf_size());
  const Eigen::VectorXd Qq_res = Eigen::VectorXd::Random(composition.Qq_size());
  const Eigen::VectorXd Qv_res = Eigen::VectorXd::Random(composition.Qv_size());
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.la() = Qa_res;
  residual.lf() = Qf_res;
  residual.lq() = Qq_res;
  residual.lv() = Qv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Fq_begin(), composition.Fq_size()).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Fv_begin(), composition.Fv_size()).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.C_begin(), composition.C_size()).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qa_begin(), composition.Qa_size()).isApprox(Qa_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qf_begin(), composition.Qf_size()).isApprox(Qf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qq_begin(), composition.Qq_size()).isApprox(Qq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(composition.Qv_begin(), composition.Qv_size()).isApprox(Qv_res));
  EXPECT_TRUE(residual.lx().head(composition.Qq_size()).isApprox(Qq_res));
  EXPECT_TRUE(residual.lx().tail(composition.Qv_size()).isApprox(Qv_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimKKT(), composition.dimKKT());
  EXPECT_EQ(residual.max_dimKKT(), composition.max_dimKKT());
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}