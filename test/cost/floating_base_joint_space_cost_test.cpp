#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/joint_space_cost.hpp"


namespace idocp {

class FloatingBaseJointSpaceCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    robot_ = Robot(urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string urdf_;
  Robot robot_;
};


TEST_F(FloatingBaseJointSpaceCostTest, zeroRefernceConstructor) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(dimq);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Zero(dimv); 
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Zero(dimv);
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(dimv);
  JointSpaceCost cost(robot_, q_weight, v_weight, a_weight, u_weight, qf_weight, 
                      vf_weight);
  Eigen::MatrixXd q_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::MatrixXd v_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv); 
  Eigen::MatrixXd a_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd u_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd qf_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::MatrixXd vf_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
  for (int i=0; i<dimq; ++i) {
    q_weight_mat(i, i) = q_weight(i);
    qf_weight_mat(i, i) = qf_weight(i);
  }
  for (int i=0; i<dimv; ++i) {
    v_weight_mat(i, i) = v_weight(i);
    a_weight_mat(i, i) = a_weight(i);
    u_weight_mat(i, i) = u_weight(i);
    vf_weight_mat(i, i) = vf_weight(i);
  }
  pinocchio::Model model;
  pinocchio::urdf::buildModel(urdf_, model);
  const Eigen::VectorXd q 
      = pinocchio::randomConfiguration(model, -Eigen::VectorXd::Ones(dimq), 
                                       Eigen::VectorXd::Ones(dimq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
  robot_.computeConfigurationJacobian(q);
  const double l_ref = 0.5 * dtau_ 
                        * ((q_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
                            + (v_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum()
                            + (a_weight.array()* (a-a_ref).array()*(a-a_ref).array()).sum()
                            + (u_weight.array()* (u-u_ref).array()*(u-u_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.l(dtau_, q, v, a, u), l_ref);
  const double phi_ref = 0.5 * ((qf_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
                                + (vf_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum());
  EXPECT_DOUBLE_EQ(cost.phi(q, v), phi_ref);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lq_ref_config = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv);
  cost.lq(robot_, dtau_, q, lq);
  lq_ref_config.array() = dtau_ * q_weight_mat * (q-q_ref);
  robot_.computeTangentGradient(lq_ref_config, lq_ref);
  EXPECT_TRUE(lq.isApprox(lq_ref));
  cost.lv(dtau_, v, lv);
  lv_ref.array() = dtau_ * v_weight_mat * (v-v_ref);
  EXPECT_TRUE(lv.isApprox(lv_ref));
  cost.la(dtau_, a, la);
  la_ref.array() = dtau_ * a_weight_mat * (a-a_ref);
  EXPECT_TRUE(la.isApprox(la_ref));
  cost.lu(dtau_, u, lu);
  lu_ref.array() = dtau_ * u_weight_mat * (u-u_ref);
  EXPECT_TRUE(lu.isApprox(lu_ref));
  cost.phiq(robot_, q, lq);
  lq_ref_config.array() = qf_weight_mat * (q-q_ref);
  robot_.computeTangentGradient(lq_ref_config, lq_ref);
  EXPECT_TRUE(lq.isApprox(lq_ref));
  std::cout << lq.transpose() << std::endl;
  std::cout << lq_ref.transpose() << std::endl;
  cost.phiv(v, lv);
  lv_ref.array() = vf_weight_mat * (v-v_ref);
  EXPECT_TRUE(lv.isApprox(lv_ref));
  Eigen::MatrixXd lqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd lqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd lvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd laa = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv, dimv);
  cost.lqq(robot_, dtau_, lqq);
  robot_.computeTangentHessian(q_weight_mat, lqq_ref);
  EXPECT_TRUE(lqq.isApprox(dtau_*lqq_ref));
  cost.lvv(dtau_, lvv);
  EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
  cost.laa(dtau_, laa);
  EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
  cost.luu(dtau_, luu);
  EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
  lqq.setZero();
  lvv.setZero();
  laa.setZero();
  luu.setZero();
  cost.augment_lqq(robot_, dtau_, lqq);
  robot_.computeTangentHessian(q_weight_mat, lqq_ref);
  EXPECT_TRUE(lqq.isApprox(dtau_*lqq_ref));
  cost.augment_lvv(dtau_, lvv);
  EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
  cost.augment_laa(dtau_, laa);
  EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
  cost.augment_luu(dtau_, luu);
  EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
  lqq.setZero();
  lvv.setZero();
  cost.phiqq(robot_, lqq);
  robot_.computeTangentHessian(qf_weight_mat, lqq_ref);
  EXPECT_TRUE(lqq.isApprox(lqq_ref));
  cost.phivv(lvv);
  EXPECT_TRUE(lvv.isApprox(vf_weight_mat));
}


// TEST_F(FloatingBaseJointSpaceCostTest, withRefernceConstructor) {
//   const int dimq = robot_.dimq();
//   const int dimv = robot_.dimv();
//   const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd q_ref = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimv);
//   JointSpaceCost cost(robot_, q_ref, v_ref, a_ref, u_ref, q_weight, v_weight, 
//                       a_weight, u_weight, qf_weight, vf_weight);
//   Eigen::MatrixXd q_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd v_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv); 
//   Eigen::MatrixXd a_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd u_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd qf_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd vf_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   for (int i=0; i<dimq; ++i) {
//     q_weight_mat(i, i) = q_weight(i);
//     qf_weight_mat(i, i) = qf_weight(i);
//   }
//   for (int i=0; i<dimv; ++i) {
//     v_weight_mat(i, i) = v_weight(i);
//     a_weight_mat(i, i) = a_weight(i);
//     u_weight_mat(i, i) = u_weight(i);
//     vf_weight_mat(i, i) = vf_weight(i);
//   }
//   pinocchio::Model model;
//   pinocchio::urdf::buildModel(urdf_, model);
//   const Eigen::VectorXd q 
//       = pinocchio::randomConfiguration(model, -Eigen::VectorXd::Ones(dimq), 
//                                        Eigen::VectorXd::Ones(dimq));
//   const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
//   robot_.computeConfigurationJacobian(q);
//   const double l_ref = 0.5 * dtau_ 
//                         * ((q_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                             + (v_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum()
//                             + (a_weight.array()* (a-a_ref).array()*(a-a_ref).array()).sum()
//                             + (u_weight.array()* (u-u_ref).array()*(u-u_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.l(dtau_, q, v, a, u), l_ref);
//   const double phi_ref = 0.5 * ((qf_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                                 + (vf_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.phi(q, v), phi_ref);
//   Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lq_ref_config = Eigen::VectorXd::Zero(dimq);
//   Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv);
//   cost.lq(robot_, dtau_, q, lq);
//   lq_ref_config.array() = dtau_ * q_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.lv(dtau_, v, lv);
//   lv_ref.array() = dtau_ * v_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   cost.la(dtau_, a, la);
//   la_ref.array() = dtau_ * a_weight_mat * (a-a_ref);
//   EXPECT_TRUE(la.isApprox(la_ref));
//   cost.lu(dtau_, u, lu);
//   lu_ref.array() = dtau_ * u_weight_mat * (u-u_ref);
//   EXPECT_TRUE(lu.isApprox(lu_ref));
//   cost.phiq(robot_, q, lq);
//   lq_ref_config.array() = qf_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.phiv(v, lv);
//   lv_ref.array() = vf_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   Eigen::MatrixXd lqq = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd lvv = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd laa = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv, dimv);
//   cost.lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   laa.setZero();
//   luu.setZero();
//   cost.augment_lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.augment_lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.augment_laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.augment_luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   cost.phiqq(robot_, lqq);
//   EXPECT_TRUE(lqq.isApprox(configuration_Jacobian.transpose()*qf_weight_mat*configuration_Jacobian));
//   cost.phivv(lvv);
//   EXPECT_TRUE(lvv.isApprox(vf_weight_mat));
// }


// TEST_F(FloatingBaseJointSpaceCostTest, setReference) {
//   const int dimq = robot_.dimq();
//   const int dimv = robot_.dimv();
//   const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd q_ref = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimv);
//   JointSpaceCost cost(robot_, q_weight, v_weight, a_weight, u_weight, qf_weight, 
//                       vf_weight);
//   cost.set_q_ref(q_ref);
//   cost.set_v_ref(v_ref);
//   cost.set_a_ref(a_ref);
//   cost.set_u_ref(u_ref);
//   Eigen::MatrixXd q_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd v_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv); 
//   Eigen::MatrixXd a_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd u_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd qf_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd vf_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   for (int i=0; i<dimq; ++i) {
//     q_weight_mat(i, i) = q_weight(i);
//     qf_weight_mat(i, i) = qf_weight(i);
//   }
//   for (int i=0; i<dimv; ++i) {
//     v_weight_mat(i, i) = v_weight(i);
//     a_weight_mat(i, i) = a_weight(i);
//     u_weight_mat(i, i) = u_weight(i);
//     vf_weight_mat(i, i) = vf_weight(i);
//   }
//   pinocchio::Model model;
//   pinocchio::urdf::buildModel(urdf_, model);
//   const Eigen::VectorXd q 
//       = pinocchio::randomConfiguration(model, -Eigen::VectorXd::Ones(dimq), 
//                                        Eigen::VectorXd::Ones(dimq));
//   const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
//   robot_.computeConfigurationJacobian(q);
//   const double l_ref = 0.5 * dtau_ 
//                         * ((q_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                             + (v_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum()
//                             + (a_weight.array()* (a-a_ref).array()*(a-a_ref).array()).sum()
//                             + (u_weight.array()* (u-u_ref).array()*(u-u_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.l(dtau_, q, v, a, u), l_ref);
//   const double phi_ref = 0.5 * ((qf_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                                 + (vf_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.phi(q, v), phi_ref);
//   Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lq_ref_config = Eigen::VectorXd::Zero(dimq);
//   Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv);
//   cost.lq(robot_, dtau_, q, lq);
//   lq_ref_config.array() = dtau_ * q_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.lv(dtau_, v, lv);
//   lv_ref.array() = dtau_ * v_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   cost.la(dtau_, a, la);
//   la_ref.array() = dtau_ * a_weight_mat * (a-a_ref);
//   EXPECT_TRUE(la.isApprox(la_ref));
//   cost.lu(dtau_, u, lu);
//   lu_ref.array() = dtau_ * u_weight_mat * (u-u_ref);
//   EXPECT_TRUE(lu.isApprox(lu_ref));
//   cost.phiq(robot_, q, lq);
//   lq_ref_config.array() = qf_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.phiv(v, lv);
//   lv_ref.array() = vf_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   Eigen::MatrixXd lqq = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd lvv = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd laa = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv, dimv);
//   cost.lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   laa.setZero();
//   luu.setZero();
//   cost.augment_lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.augment_lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.augment_laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.augment_luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   cost.phiqq(robot_, lqq);
//   EXPECT_TRUE(lqq.isApprox(configuration_Jacobian.transpose()*qf_weight_mat*configuration_Jacobian));
//   cost.phivv(lvv);
//   EXPECT_TRUE(lvv.isApprox(vf_weight_mat));
// }


// TEST_F(FloatingBaseJointSpaceCostTest, setWeights) {
//   const int dimq = robot_.dimq();
//   const int dimv = robot_.dimv();
//   const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd q_ref = Eigen::VectorXd::Random(dimq);
//   const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
//   const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimv);
//   JointSpaceCost cost(robot_, Eigen::VectorXd::Zero(dimq), 
//                       Eigen::VectorXd::Zero(dimv), Eigen::VectorXd::Zero(dimv), 
//                       Eigen::VectorXd::Zero(dimv), Eigen::VectorXd::Zero(dimq), 
//                       Eigen::VectorXd::Zero(dimv));
//   cost.set_q_weight(q_weight);
//   cost.set_v_weight(v_weight);
//   cost.set_a_weight(a_weight);
//   cost.set_u_weight(u_weight);
//   cost.set_qf_weight(qf_weight);
//   cost.set_vf_weight(vf_weight);
//   cost.set_q_ref(q_ref);
//   cost.set_v_ref(v_ref);
//   cost.set_a_ref(a_ref);
//   cost.set_u_ref(u_ref);
//   Eigen::MatrixXd q_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd v_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv); 
//   Eigen::MatrixXd a_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd u_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd qf_weight_mat = Eigen::MatrixXd::Zero(dimq, dimq);
//   Eigen::MatrixXd vf_weight_mat = Eigen::MatrixXd::Zero(dimv, dimv);
//   for (int i=0; i<dimq; ++i) {
//     q_weight_mat(i, i) = q_weight(i);
//     qf_weight_mat(i, i) = qf_weight(i);
//   }
//   for (int i=0; i<dimv; ++i) {
//     v_weight_mat(i, i) = v_weight(i);
//     a_weight_mat(i, i) = a_weight(i);
//     u_weight_mat(i, i) = u_weight(i);
//     vf_weight_mat(i, i) = vf_weight(i);
//   }
//   pinocchio::Model model;
//   pinocchio::urdf::buildModel(urdf_, model);
//   const Eigen::VectorXd q 
//       = pinocchio::randomConfiguration(model, -Eigen::VectorXd::Ones(dimq), 
//                                        Eigen::VectorXd::Ones(dimq));
//   const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
//   const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
//   robot_.computeConfigurationJacobian(q);
//   const double l_ref = 0.5 * dtau_ 
//                         * ((q_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                             + (v_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum()
//                             + (a_weight.array()* (a-a_ref).array()*(a-a_ref).array()).sum()
//                             + (u_weight.array()* (u-u_ref).array()*(u-u_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.l(dtau_, q, v, a, u), l_ref);
//   const double phi_ref = 0.5 * ((qf_weight.array()* (q-q_ref).array()*(q-q_ref).array()).sum()
//                                 + (vf_weight.array()* (v-v_ref).array()*(v-v_ref).array()).sum());
//   EXPECT_DOUBLE_EQ(cost.phi(q, v), phi_ref);
//   Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lq_ref_config = Eigen::VectorXd::Zero(dimq);
//   Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lv_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd la_ref = Eigen::VectorXd::Zero(dimv);
//   Eigen::VectorXd lu_ref = Eigen::VectorXd::Zero(dimv);
//   cost.lq(robot_, dtau_, q, lq);
//   lq_ref_config.array() = dtau_ * q_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.lv(dtau_, v, lv);
//   lv_ref.array() = dtau_ * v_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   cost.la(dtau_, a, la);
//   la_ref.array() = dtau_ * a_weight_mat * (a-a_ref);
//   EXPECT_TRUE(la.isApprox(la_ref));
//   cost.lu(dtau_, u, lu);
//   lu_ref.array() = dtau_ * u_weight_mat * (u-u_ref);
//   EXPECT_TRUE(lu.isApprox(lu_ref));
//   cost.phiq(robot_, q, lq);
//   lq_ref_config.array() = qf_weight_mat * (q-q_ref);
//   lq_ref = configuration_Jacobian.transpose() * lq_ref_config;
//   EXPECT_TRUE(lq.isApprox(lq_ref));
//   cost.phiv(v, lv);
//   lv_ref.array() = vf_weight_mat * (v-v_ref);
//   EXPECT_TRUE(lv.isApprox(lv_ref));
//   Eigen::MatrixXd lqq = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd lvv = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd laa = Eigen::MatrixXd::Zero(dimv, dimv);
//   Eigen::MatrixXd luu = Eigen::MatrixXd::Zero(dimv, dimv);
//   cost.lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   laa.setZero();
//   luu.setZero();
//   cost.augment_lqq(robot_, dtau_, lqq);
//   EXPECT_TRUE(lqq.isApprox(dtau_*configuration_Jacobian.transpose()*q_weight_mat*configuration_Jacobian));
//   cost.augment_lvv(dtau_, lvv);
//   EXPECT_TRUE(lvv.isApprox(dtau_*v_weight_mat));
//   cost.augment_laa(dtau_, laa);
//   EXPECT_TRUE(laa.isApprox(dtau_*a_weight_mat));
//   cost.augment_luu(dtau_, luu);
//   EXPECT_TRUE(luu.isApprox(dtau_*u_weight_mat));
//   lqq.setZero();
//   lvv.setZero();
//   cost.phiqq(robot_, lqq);
//   EXPECT_TRUE(lqq.isApprox(configuration_Jacobian.transpose()*qf_weight_mat*configuration_Jacobian));
//   cost.phivv(lvv);
//   EXPECT_TRUE(lvv.isApprox(vf_weight_mat));
// }


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}