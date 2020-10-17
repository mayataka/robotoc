#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"


namespace idocp {

class ImpulseKKTResidualTest : public ::testing::Test {
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


TEST_F(ImpulseKKTResidualTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus({true});
  ImpulseKKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = 2*contact_status.dimf();
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimf+dimc);
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  const Eigen::VectorXd lf_res = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd lq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd ldv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv_res = Eigen::VectorXd::Random(dimv);
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.lf() = lf_res;
  residual.lq() = lq_res;
  residual.lv() = lv_res;
  residual.ldv = ldv_res;
  residual.dv_res = dv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(               0, dimv).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(            dimv, dimv).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(          2*dimv, dimc).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     2*dimv+dimc, dimf).isApprox(lf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(2*dimv+dimc+dimf, dimv).isApprox(lq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(3*dimv+dimc+dimf, dimv).isApprox(lv_res));
  EXPECT_TRUE(residual.C_contact_position().isApprox(C_res.head(dimf)));
  EXPECT_TRUE(residual.C_contact_velocity().isApprox(C_res.tail(dimf)));
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_TRUE(residual.lx().head(dimv).isApprox(lq_res));
  EXPECT_TRUE(residual.lx().tail(dimv).isApprox(lv_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimc+dimf);
  EXPECT_EQ(residual.max_dimKKT(), 4*dimv+3*robot.max_dimf());
}


TEST_F(ImpulseKKTResidualTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  ImpulseKKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = 2*contact_status.dimf();
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimf+dimc);
  const Eigen::VectorXd Fq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd Fv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd C_res = Eigen::VectorXd::Random(dimc);
  const Eigen::VectorXd lf_res = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd lq_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd ldv_res = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv_res = Eigen::VectorXd::Random(dimv);
  residual.Fq() = Fq_res;
  residual.Fv() = Fv_res;
  residual.C() = C_res;
  residual.lf() = lf_res;
  residual.lq() = lq_res;
  residual.lv() = lv_res;
  residual.ldv = ldv_res;
  residual.dv_res = dv_res;
  EXPECT_TRUE(residual.KKT_residual().segment(               0, dimv).isApprox(Fq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(            dimv, dimv).isApprox(Fv_res));
  EXPECT_TRUE(residual.KKT_residual().segment(          2*dimv, dimc).isApprox(C_res));
  EXPECT_TRUE(residual.KKT_residual().segment(     2*dimv+dimc, dimf).isApprox(lf_res));
  EXPECT_TRUE(residual.KKT_residual().segment(2*dimv+dimc+dimf, dimv).isApprox(lq_res));
  EXPECT_TRUE(residual.KKT_residual().segment(3*dimv+dimc+dimf, dimv).isApprox(lv_res));
  EXPECT_TRUE(residual.C_contact_position().isApprox(C_res.head(dimf)));
  EXPECT_TRUE(residual.C_contact_velocity().isApprox(C_res.tail(dimf)));
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_TRUE(residual.lx().head(dimv).isApprox(lq_res));
  EXPECT_TRUE(residual.lx().tail(dimv).isApprox(lv_res));
  residual.setZero();
  EXPECT_TRUE(residual.KKT_residual().isZero());
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimc(), dimc);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimc+dimf);
  EXPECT_EQ(residual.max_dimKKT(), 4*dimv+3*robot.max_dimf());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}