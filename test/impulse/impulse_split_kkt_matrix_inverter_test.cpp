#include <random>
#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hpp"


namespace idocp {

class ImpulseSplitKKTMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  std::string fixed_base_urdf, floating_base_urdf;
};


TEST_F(ImpulseSplitKKTMatrixInverterTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  ImpulseStatus impulse_status(contact_frames.size());
  impulse_status.setImpulseStatus({true});
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimf();
  const int dimQ = 2*dimv + dimf;
  const int dimKKT = 4*dimv + 2*dimf;
  const Eigen::MatrixXd KKT_seed_mat = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd KKT_mat = KKT_seed_mat * KKT_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  KKT_mat.topLeftCorner(dimx+dimf, dimx+dimf).setZero();
  KKT_mat.block(               0,        dimx+dimf, dimv, dimf).setZero(); // Fqf
  KKT_mat.block(               0,      dimx+2*dimf, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv); // Fqq
  KKT_mat.block(               0, dimx+2*dimf+dimv, dimv, dimv).setZero(); // Fqv
  KKT_mat.block(            dimv, dimx+2*dimf+dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv); // Fvv
  KKT_mat.block(            dimx,        dimx+dimf, dimf, dimf).setZero(); // Vf
  KKT_mat.template triangularView<Eigen::StrictlyLower>() 
      = KKT_mat.transpose().template triangularView<Eigen::StrictlyLower>();
  ImpulseSplitKKTMatrixInverter inverter(robot);
  const Eigen::MatrixXd FC = KKT_mat.topRightCorner(dimQ, dimQ);
  const Eigen::MatrixXd multiplied_mat = Eigen::MatrixXd::Random(dimQ, dimQ);
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(dimQ, dimQ);
  inverter.multiplyJac(FC, multiplied_mat, res);
  const Eigen::MatrixXd res_ref = FC * multiplied_mat;
  EXPECT_TRUE(res.isApprox(res_ref));
  const Eigen::MatrixXd Q = KKT_mat.bottomRightCorner(dimQ, dimQ);
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  inverter.invert(FC, Q, KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref, 1.0e-08));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity(1.0e-06));
}


TEST_F(ImpulseSplitKKTMatrixInverterTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  auto impulse_status = robot.createImpulseStatus();
  std::random_device rnd;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (rnd() % 2 == 0) {
      impulse_status.activateImpulse(i);
    }
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimf();
  const int dimQ = 2*dimv + dimf;
  const int dimKKT = 4*dimv + 2*dimf;
  const Eigen::MatrixXd KKT_seed_mat = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd KKT_mat = KKT_seed_mat * KKT_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  KKT_mat.topLeftCorner(dimx+dimf, dimx+dimf).setZero();
  KKT_mat.block(               0,        dimx+dimf, dimv, dimf).setZero(); // Fqf
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Zero(dimv, dimv);
  robot.dSubtractdConfigurationMinus(robot.generateFeasibleConfiguration(), 
                                     robot.generateFeasibleConfiguration(), Fqq);
  KKT_mat.block(               0,      dimx+2*dimf, dimv, dimv) = Fqq; // Fqq
  KKT_mat.block(               0, dimx+2*dimf+dimv, dimv, dimv).setZero(); // Fqv
  KKT_mat.block(            dimv, dimx+2*dimf+dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv); // Fvv
  KKT_mat.block(            dimx,        dimx+dimf, dimf, dimf).setZero(); // Vf
  KKT_mat.template triangularView<Eigen::StrictlyLower>() 
      = KKT_mat.transpose().template triangularView<Eigen::StrictlyLower>();
  ImpulseSplitKKTMatrixInverter inverter(robot);
  const Eigen::MatrixXd FC = KKT_mat.topRightCorner(dimQ, dimQ);
  const Eigen::MatrixXd multiplied_mat = Eigen::MatrixXd::Random(dimQ, dimQ);
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(dimQ, dimQ);
  inverter.multiplyJac(FC, multiplied_mat, res);
  const Eigen::MatrixXd res_ref = FC * multiplied_mat;
  EXPECT_TRUE(res.isApprox(res_ref));
  const Eigen::MatrixXd Q = KKT_mat.bottomRightCorner(dimQ, dimQ);
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  inverter.invert(FC, Q, KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref, 1.0e-08));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity(1.0e-06));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}