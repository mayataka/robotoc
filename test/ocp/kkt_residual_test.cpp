#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_residual.hpp"


namespace idocp {

class KKTResidualTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
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
  robot.setContactStatus(contact_status);
  KKTResidual residual(robot);
  residual.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  const Eigen::VectorXd Qa_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Qf_res = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd Qq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Qv_res = Eigen::VectorXd::Random(dimv);
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.la() = Qa_res;
  residual.lf() = Qf_res;
  residual.lq() = Qq_res;
  residual.lv() = Qv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(               0, dimv).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(            dimv, dimv).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(          2*dimv, dimc).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     2*dimv+dimc, dimv).isApprox(Qa_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     3*dimv+dimc, dimf).isApprox(Qf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(3*dimv+dimc+dimf, dimv).isApprox(Qq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(4*dimv+dimc+dimf, dimv).isApprox(Qv_res));
  EXPECT_TRUE(residual.C_floating_base().isApprox(C_res.head(robot.dim_passive())));
  EXPECT_TRUE(residual.C_contacts().isApprox(C_res.tail(robot.dimf())));
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_TRUE(residual.lx().head(dimv).isApprox(Qq_res));
  EXPECT_TRUE(residual.lx().tail(dimv).isApprox(Qv_res));
  EXPECT_EQ(residual.laf().size(), dimv+dimf);
  EXPECT_TRUE(residual.laf().head(dimv).isApprox(Qa_res));
  EXPECT_TRUE(residual.laf().tail(dimf).isApprox(Qf_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 5*dimv+dimc+dimf);
  EXPECT_EQ(residual.max_dimKKT(), 5*dimv+2*robot.max_dimf()+robot.dim_passive());
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
  KKTResidual residual(robot);
  residual.setContactStatus(robot);
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dim_passive = robot.dim_passive();
  const int dimc = robot.dim_passive() + robot.dimf();
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  const Eigen::VectorXd Qa_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Qf_res = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd Qq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Qv_res = Eigen::VectorXd::Random(dimv);
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.la() = Qa_res;
  residual.lf() = Qf_res;
  residual.lq() = Qq_res;
  residual.lv() = Qv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(               0, dimv).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(            dimv, dimv).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(          2*dimv, dimc).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     2*dimv+dimc, dimv).isApprox(Qa_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     3*dimv+dimc, dimf).isApprox(Qf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(3*dimv+dimc+dimf, dimv).isApprox(Qq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(4*dimv+dimc+dimf, dimv).isApprox(Qv_res));
  EXPECT_TRUE(residual.C_floating_base().isApprox(C_res.head(robot.dim_passive())));
  EXPECT_TRUE(residual.C_contacts().isApprox(C_res.tail(robot.dimf())));
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_TRUE(residual.lx().head(dimv).isApprox(Qq_res));
  EXPECT_TRUE(residual.lx().tail(dimv).isApprox(Qv_res));
  EXPECT_EQ(residual.laf().size(), dimv+dimf);
  EXPECT_TRUE(residual.laf().head(dimv).isApprox(Qa_res));
  EXPECT_TRUE(residual.laf().tail(dimf).isApprox(Qf_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 5*dimv+dimc+dimf);
  EXPECT_EQ(residual.max_dimKKT(), 5*dimv+2*robot.max_dimf()+robot.dim_passive());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}