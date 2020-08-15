#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/pdipm_func.hpp"


namespace idocp {

class JointPositionLowerLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    barrier_ = 1.0e-08;
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double barrier_, dtau_;
  Eigen::VectorXd slack_, dual_, dslack_, ddual_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(JointPositionLowerLimitTest, isFeasibleFixedBase) {
  JointPositionLowerLimit limit(fixed_base_robot_); 
  ConstraintComponentData data(limit.dimc());
  Eigen::VectorXd q(fixed_base_robot_.dimq());
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(fixed_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  fixed_base_robot_.generateFeasibleConfiguration(q);
  EXPECT_TRUE(limit.isFeasible(fixed_base_robot_, data, a, f, q, v, u));
  q = 2*fixed_base_robot_.lowerJointPositionLimit();
  EXPECT_FALSE(limit.isFeasible(fixed_base_robot_, data, a, f, q, v, u));
}


TEST_F(JointPositionLowerLimitTest, isFeasibleFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  Eigen::VectorXd q(floating_base_robot_.dimq());
  Eigen::VectorXd v = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(floating_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  floating_base_robot_.generateFeasibleConfiguration(q);
  EXPECT_TRUE(limit.isFeasible(floating_base_robot_, data, a, f, q, v, u));
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  q.tail(dimc) = 2*floating_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(q.size(), floating_base_robot_.dimq());
  EXPECT_FALSE(limit.isFeasible(floating_base_robot_, data, a, f, q, v, u));
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFixedBase) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  Eigen::VectorXd q(dimq);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(fixed_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(q);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, la, lf, lq, lv);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, lu);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (q-qmin);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  pdipm::pdipmfunc::SetSlackAndDualPositive(dimq, barrier_, slack_ref, dual_ref);
  lq_ref = - dtau_ * dual_ref;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimq, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  const double l1residual_ref = (dtau_*(qmin-q)+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.residualSquaredNrom(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipm::pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(qmin-q)+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  Eigen::VectorXd q(dimq);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(fixed_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd qmin = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  ASSERT_EQ(dimc+6, dimv);
  floating_base_robot_.generateFeasibleConfiguration(q);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, la, lf, lq, lv);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, lu);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (q.tail(dimc)-qmin);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  pdipm::pdipmfunc::SetSlackAndDualPositive(dimc, barrier_, slack_ref, dual_ref);
  lq_ref.tail(dimc) = - dtau_ * dual_ref;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimc, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  const double l1residual_ref = (dtau_*(qmin-q.tail(dimc))+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.residualSquaredNrom(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipm::pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(qmin-q.tail(dimc))+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFixedBase) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  Eigen::VectorXd q(dimq);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(fixed_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(q);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::MatrixXd Cvv = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::MatrixXd Caa = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::MatrixXd Cff = Eigen::MatrixXd::Zero(dimf, dimf);
  Eigen::MatrixXd Cuu = Eigen::MatrixXd::Zero(dimq, dimq);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (q-qmin);
  pdipm::pdipmfunc::SetSlackAndDualPositive(dimq, barrier_, slack_ref, dual_ref);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, a, f, q, v, 
                             Caa, Cff, Cqq, Cvv, la, lf, lq, lv);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, u, Cuu, lu);
  Eigen::VectorXd residual_ref = dtau_ * (qmin-q) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipm::pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  lq_ref.array() 
      -= dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Zero(dimq, dimq);
  for (int i=0; i<dimq; ++i) {
    Cqq_ref(i, i) += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lu.isZero());
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isZero());
  EXPECT_TRUE(Caa.isZero());
  EXPECT_TRUE(Cff.isZero());
  EXPECT_TRUE(Cuu.isZero());
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd da = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd df = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(fixed_base_robot_, data, dtau_,  
                                     da, df, dq, dv, du);
  const Eigen::VectorXd dslack_ref = dtau_ * dq - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimq);
  pdipm::pdipmfunc::ComputeDualDirection(dual_ref, slack_ref, dslack_ref,  
                                         duality_ref, ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref 
      = pdipm::pdipmfunc::FractionToBoundary(dimq, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipm::pdipmfunc::FractionToBoundary(dimq, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimq, barrier_, 
                                           slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  const double cost_slack_barrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimq, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  lq.setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, la, lf, lq, lv);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, lu);
  lq_ref = - dtau_ * dual_ref;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(lu.isZero());
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  Eigen::VectorXd q(dimq);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  Eigen::VectorXd f = Eigen::VectorXd::Zero(floating_base_robot_.dimf());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  Eigen::VectorXd qmin = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  floating_base_robot_.generateFeasibleConfiguration(q);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, a, f, q, v, u);
  Eigen::VectorXd la = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  Eigen::VectorXd lq = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lv = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd lu = Eigen::VectorXd::Zero(dimv);
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Cvv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Caa = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Cff = Eigen::MatrixXd::Zero(dimf, dimf);
  Eigen::MatrixXd Cuu = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (q.tail(dimc)-qmin);
  pdipm::pdipmfunc::SetSlackAndDualPositive(dimc, barrier_, slack_ref, dual_ref);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, a, f, q, v, 
                             Caa, Cff, Cqq, Cvv, la, lf, lq, lv);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, u, Cuu, lu);
  Eigen::VectorXd residual_ref = dtau_ * (qmin-q.tail(dimc)) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipm::pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  lq_ref.tail(dimc).array() 
      -= dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  for (int i=0; i<dimc; ++i) {
    Cqq_ref(6+i, 6+i) += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lu.isZero());
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isZero());
  EXPECT_TRUE(Caa.isZero());
  EXPECT_TRUE(Cff.isZero());
  EXPECT_TRUE(Cuu.isZero());
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd da = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd df = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(floating_base_robot_, data, dtau_,  
                                     da, df, dq, dv, du);
  const Eigen::VectorXd dslack_ref = dtau_ * dq.tail(dimc) - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimc);
  pdipm::pdipmfunc::ComputeDualDirection(dual_ref, slack_ref, dslack_ref,  
                                         duality_ref, ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref 
      = pdipm::pdipmfunc::FractionToBoundary(dimc, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipm::pdipmfunc::FractionToBoundary(dimc, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimc, barrier_, 
                                           slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  const double cost_slack_barrier_ref 
      = pdipm::pdipmfunc::SlackBarrierCost(dimc, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  lq.setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, la, lf, lq, lv);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, lu);
  lq_ref.tail(dimc) = - dtau_ * dual_ref;
  EXPECT_TRUE(lq.isApprox(lq_ref));
  EXPECT_TRUE(la.isZero());
  EXPECT_TRUE(lf.isZero());
  EXPECT_TRUE(lv.isZero());
  EXPECT_TRUE(lu.isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}