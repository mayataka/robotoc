#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/pdipm_func.hpp"

#include "joint_variables_lower_limits.hpp"

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
  SplitSolution s(fixed_base_robot_);
  fixed_base_robot_.normalizeConfiguration(s.q);
  EXPECT_TRUE(limit.isFeasible(fixed_base_robot_, data, s));
  s.q = 2*fixed_base_robot_.lowerJointPositionLimit();
  EXPECT_FALSE(limit.isFeasible(fixed_base_robot_, data, s));
}


TEST_F(JointPositionLowerLimitTest, isFeasibleFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  SplitSolution s(floating_base_robot_);
  floating_base_robot_.normalizeConfiguration(s.q);
  EXPECT_TRUE(limit.isFeasible(floating_base_robot_, data, s));
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  s.q.tail(dimc) = 2*floating_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(s.q.size(), floating_base_robot_.dimq());
  EXPECT_FALSE(limit.isFeasible(floating_base_robot_, data, s));
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFixedBase) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  SplitSolution s(fixed_base_robot_);
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(s.q);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual.lu);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (s.q-qmin);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack_ref, dual_ref);
  lq_ref = - dtau_ * dual_ref;
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref = pdipmfunc::CostSlackBarrier(barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(fixed_base_robot_, data, dtau_, s);
  const double l1residual_ref = (dtau_*(qmin-s.q)+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.squaredKKTErrorNorm(fixed_base_robot_, data, dtau_, s);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDuality(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(qmin-s.q)+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFixedBase2) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  SplitSolution s(fixed_base_robot_);
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(s.q);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual.lu);
  pdipmtest::JointVariablesLowerLimits limit_ref(fixed_base_robot_, fixed_base_robot_.lowerJointPositionLimit(), barrier_);
  limit_ref.setSlackAndDual(dtau_, s.q);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  limit_ref.augmentDualResidual(dtau_, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref = limit_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(fixed_base_robot_, data, dtau_, s);
  const double l1residual_ref = limit_ref.residualL1Nrom(dtau_, s.q);
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.squaredKKTErrorNorm(fixed_base_robot_, data, dtau_, s);
  const double l2residual_ref = limit_ref.residualSquaredNrom(dtau_, s.q);
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  SplitSolution s(floating_base_robot_);
  Eigen::VectorXd qmin = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  ASSERT_EQ(dimc+6, dimv);
  floating_base_robot_.generateFeasibleConfiguration(s.q);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, s);
  KKTMatrix kkt_matrix(floating_base_robot_);
  KKTResidual kkt_residual(floating_base_robot_);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual.lu);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (s.q.tail(dimc)-qmin);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack_ref, dual_ref);
  lq_ref.tail(dimc) = - dtau_ * dual_ref;
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref 
      = pdipmfunc::CostSlackBarrier(barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(floating_base_robot_, data, dtau_, s);
  const double l1residual_ref = (dtau_*(qmin-s.q.tail(dimc))+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.squaredKKTErrorNorm(floating_base_robot_, data, dtau_, s);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDuality(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(qmin-s.q.tail(dimc))+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDualFloatingBase2) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  SplitSolution s(floating_base_robot_);
  Eigen::VectorXd qmin = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  ASSERT_EQ(dimc+6, dimv);
  floating_base_robot_.generateFeasibleConfiguration(s.q);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, s);
  KKTMatrix kkt_matrix(floating_base_robot_);
  KKTResidual kkt_residual(floating_base_robot_);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual.lu);
  pdipmtest::JointVariablesLowerLimits limit_ref(floating_base_robot_, floating_base_robot_.lowerJointPositionLimit(), barrier_);
  limit_ref.setSlackAndDual(dtau_, s.q);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  limit_ref.augmentDualResidual(dtau_, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref = limit_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(fixed_base_robot_, data, dtau_, s);
  const double l1residual_ref = limit_ref.residualL1Nrom(dtau_, s.q);
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.squaredKKTErrorNorm(fixed_base_robot_, data, dtau_, s);
  const double l2residual_ref = limit_ref.residualSquaredNrom(dtau_, s.q);
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFixedBase) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  SplitSolution s(fixed_base_robot_);
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(s.q);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (s.q-qmin);
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack_ref, dual_ref);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  Eigen::VectorXd residual_ref = dtau_ * (qmin-s.q) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDuality(barrier_, slack_ref, dual_ref, duality_ref);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  lq_ref.array() 
      -= dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimq, dimq);
  for (int i=0; i<dimq; ++i) {
    Qqq_ref(i, i) += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(Qqq_ref));
  EXPECT_TRUE(kkt_matrix.Qvv().isZero());
  EXPECT_TRUE(kkt_matrix.Qaa().isZero());
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Quu.isZero());
  SplitDirection d(fixed_base_robot_);
  d.dq() = Eigen::VectorXd::Random(dimv);
  d.dv() = Eigen::VectorXd::Random(dimv);
  d.da() = Eigen::VectorXd::Random(dimv);
  d.df() = Eigen::VectorXd::Random(dimf);
  d.du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(fixed_base_robot_, data, dtau_, d);
  const Eigen::VectorXd dslack_ref = dtau_ * d.dq() - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDualDirection(slack_ref, dual_ref, dslack_ref, duality_ref, 
                                  ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimq, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimq, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref 
      = pdipmfunc::CostSlackBarrier(barrier_, slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  const double cost_slack_barrier_ref 
      = pdipmfunc::CostSlackBarrier(barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  kkt_residual.lq().setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual.lu);
  lq_ref = - dtau_ * dual_ref;
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFixedBase2) {
  JointPositionLowerLimit limit(fixed_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = fixed_base_robot_.dimq();
  const int dimv = fixed_base_robot_.dimv();
  const int dimf = fixed_base_robot_.dimf();
  SplitSolution s(fixed_base_robot_);
  Eigen::VectorXd qmin = fixed_base_robot_.lowerJointPositionLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.lowerJointPositionLimit().size());
  fixed_base_robot_.generateFeasibleConfiguration(s.q);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  limit.setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  pdipmtest::JointVariablesLowerLimits limit_ref(fixed_base_robot_, fixed_base_robot_.lowerJointPositionLimit(), barrier_);
  limit_ref.setSlackAndDual(dtau_, s.q);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  limit.condenseSlackAndDual(fixed_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimq, dimq);
  limit_ref.condenseSlackAndDual(dtau_, s.q, Qqq_ref, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(Qqq_ref));
  EXPECT_TRUE(kkt_matrix.Qvv().isZero());
  EXPECT_TRUE(kkt_matrix.Qaa().isZero());
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Quu.isZero());
  SplitDirection d(fixed_base_robot_);
  d.dq() = Eigen::VectorXd::Random(dimv);
  d.dv() = Eigen::VectorXd::Random(dimv);
  d.da() = Eigen::VectorXd::Random(dimv);
  d.df() = Eigen::VectorXd::Random(dimf);
  d.du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(fixed_base_robot_, data, dtau_, d);
  limit_ref.computeSlackAndDualDirection(dtau_, d.dq());
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref = limit_ref.maxSlackStepSize(margin_rate);
  const double dual_step_size_ref = limit_ref.maxDualStepSize(margin_rate);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref = limit_ref.costSlackBarrier(step_size);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  limit_ref.updateSlack(step_size);
  limit_ref.updateDual(step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref = limit_ref.costSlackBarrier(); 
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  kkt_residual.lq().setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(fixed_base_robot_, data, dtau_, kkt_residual.lu);
  limit_ref.augmentDualResidual(dtau_, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFloatingBase) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  SplitSolution s(floating_base_robot_);
  Eigen::VectorXd qmin = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  floating_base_robot_.generateFeasibleConfiguration(s.q);
  KKTMatrix kkt_matrix(floating_base_robot_);
  KKTResidual kkt_residual(floating_base_robot_);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, s);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (s.q.tail(dimc)-qmin);
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack_ref, dual_ref);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  Eigen::VectorXd residual_ref = dtau_ * (qmin-s.q.tail(dimc)) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDuality(barrier_, slack_ref, dual_ref, duality_ref);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  lq_ref.tail(dimc).array() 
      -= dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  for (int i=0; i<dimc; ++i) {
    Qqq_ref(6+i, 6+i) += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(Qqq_ref));
  EXPECT_TRUE(kkt_matrix.Qvv().isZero());
  EXPECT_TRUE(kkt_matrix.Qaa().isZero());
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Quu.isZero());
  SplitDirection d(floating_base_robot_);
  d.dq() = Eigen::VectorXd::Random(dimv);
  d.dv() = Eigen::VectorXd::Random(dimv);
  d.da() = Eigen::VectorXd::Random(dimv);
  d.df() = Eigen::VectorXd::Random(dimf);
  d.du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(floating_base_robot_, data, dtau_, d);
  const Eigen::VectorXd dslack_ref = dtau_ * d.dq().tail(dimc) - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDualDirection(slack_ref, dual_ref, dslack_ref, duality_ref, 
                                  ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimc, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimc, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref 
      = pdipmfunc::CostSlackBarrier(barrier_, slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  const double cost_slack_barrier_ref 
      = pdipmfunc::CostSlackBarrier(barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  kkt_residual.lq().setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual.lu);
  lq_ref.tail(dimc) = - dtau_ * dual_ref;
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDualFloatingBase2) {
  JointPositionLowerLimit limit(floating_base_robot_);
  ConstraintComponentData data(limit.dimc());
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dimf = floating_base_robot_.dimf();
  SplitSolution s(floating_base_robot_);
  Eigen::VectorXd qmax = floating_base_robot_.lowerJointPositionLimit();
  const int dimc = floating_base_robot_.lowerJointPositionLimit().size();
  floating_base_robot_.generateFeasibleConfiguration(s.q);
  KKTMatrix kkt_matrix(floating_base_robot_);
  KKTResidual kkt_residual(floating_base_robot_);
  limit.setSlackAndDual(floating_base_robot_, data, dtau_, s);
  pdipmtest::JointVariablesLowerLimits limit_ref(floating_base_robot_, floating_base_robot_.lowerJointPositionLimit(), barrier_);
  limit_ref.setSlackAndDual(dtau_, s.q);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  limit.condenseSlackAndDual(floating_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  Eigen::VectorXd lq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::MatrixXd Qqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  limit_ref.condenseSlackAndDual(dtau_, s.q, Qqq_ref, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
  EXPECT_TRUE(kkt_matrix.Qqq().isApprox(Qqq_ref));
  EXPECT_TRUE(kkt_matrix.Qvv().isZero());
  EXPECT_TRUE(kkt_matrix.Qaa().isZero());
  EXPECT_TRUE(kkt_matrix.Qff().isZero());
  EXPECT_TRUE(kkt_matrix.Quu.isZero());
  SplitDirection d(floating_base_robot_);
  d.dq() = Eigen::VectorXd::Random(dimv);
  d.dv() = Eigen::VectorXd::Random(dimv);
  d.da() = Eigen::VectorXd::Random(dimv);
  d.df() = Eigen::VectorXd::Random(dimf);
  d.du = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(floating_base_robot_, data, dtau_, d);
  limit_ref.computeSlackAndDualDirection(dtau_, d.dq());
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(data);
  const double dual_step_size = limit.maxDualStepSize(data);
  const double slack_step_size_ref = limit_ref.maxSlackStepSize(margin_rate);
  const double dual_step_size_ref = limit_ref.maxDualStepSize(margin_rate);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(data, step_size);
  const double berrier_ref = limit_ref.costSlackBarrier(step_size);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(data, step_size);
  limit.updateDual(data, step_size);
  limit_ref.updateSlack(step_size);
  limit_ref.updateDual(step_size);
  const double cost_slack_barrier = limit.costSlackBarrier(data);
  const double cost_slack_barrier_ref = limit_ref.costSlackBarrier(); 
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  kkt_residual.lq().setZero();
  lq_ref.setZero();
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual);
  limit.augmentDualResidual(floating_base_robot_, data, dtau_, kkt_residual.lu);
  limit_ref.augmentDualResidual(dtau_, lq_ref);
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_TRUE(kkt_residual.la().isZero());
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu.isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}