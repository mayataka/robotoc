#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/pdipm_func.hpp"
#include "constraints/joint_space_constraints/joint_variables_upper_limits.hpp"


namespace idocp {
namespace pdipm {

class JointVelocityUpperLimitsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf_ = "../../../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../../../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    barrier_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    while (barrier_ == 0) {
      barrier_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    }
  }

  virtual void TearDown() {
  }

  double barrier_, dtau_;
  Eigen::VectorXd slack_, dual_, dslack_, ddual_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(JointVelocityUpperLimitsTest, isFeasibleFixedBase) {
  JointVariablesUpperLimits limit(fixed_base_robot_, 
                                  fixed_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  EXPECT_TRUE(limit.isFeasible(v));
}


TEST_F(JointVelocityUpperLimitsTest, isFeasibleFloatingBase) {
  JointVariablesUpperLimits limit(floating_base_robot_, 
                                  floating_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(floating_base_robot_.dimv());
  EXPECT_TRUE(limit.isFeasible(v));
}


TEST_F(JointVelocityUpperLimitsTest, setSlackAndDualFixedBase) {
  JointVariablesUpperLimits limit(fixed_base_robot_, 
                                  fixed_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  const int dimq = fixed_base_robot_.dimq();
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd vmax = fixed_base_robot_.jointVelocityLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.jointVelocityLimit().size());
  limit.setSlackAndDual(dtau_, v);
  Eigen::VectorXd Cv = Eigen::VectorXd::Zero(dimq);
  limit.augmentDualResidual(dtau_, Cv);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (vmax-v);
  Eigen::VectorXd Cv_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::SetSlackAndDualPositive(dimq, barrier_, slack_ref, dual_ref);
  Cv_ref = dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  const double cost_slack_barrier = limit.costSlackBarrier();
  const double cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimq, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(dtau_, v);
  const double l1residual_ref = (dtau_*(v-vmax)+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.residualSquaredNrom(dtau_, v);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(v-vmax)+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointVelocityUpperLimitsTest, setSlackAndDualFloatingBase) {
  JointVariablesUpperLimits limit(floating_base_robot_, 
                                  floating_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  Eigen::VectorXd v = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd vmax = floating_base_robot_.jointVelocityLimit();
  const int dimc = floating_base_robot_.jointVelocityLimit().size();
  ASSERT_EQ(dimc+6, dimv);
  limit.setSlackAndDual(dtau_, v);
  Eigen::VectorXd Cv = Eigen::VectorXd::Zero(dimv);
  limit.augmentDualResidual(dtau_, Cv);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (vmax-v.tail(dimc));
  Eigen::VectorXd Cv_ref = Eigen::VectorXd::Zero(dimv);
  pdipmfunc::SetSlackAndDualPositive(dimc, barrier_, slack_ref, dual_ref);
  Cv_ref.tail(dimc) = dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  const double cost_slack_barrier = limit.costSlackBarrier();
  const double cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimc, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double l1residual = limit.residualL1Nrom(dtau_, v);
  const double l1residual_ref = (dtau_*(v.tail(dimc)-vmax)+slack_ref).lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1residual, l1residual_ref);
  const double l2residual = limit.residualSquaredNrom(dtau_, v);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  const double l2residual_ref 
      = (dtau_*(v.tail(dimc)-vmax)+slack_ref).squaredNorm() + duality_ref.squaredNorm();
  EXPECT_DOUBLE_EQ(l2residual, l2residual_ref);
}


TEST_F(JointVelocityUpperLimitsTest, condenseSlackAndDualFixedBase) {
  JointVariablesUpperLimits limit(fixed_base_robot_, 
                                  fixed_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  const int dimq = fixed_base_robot_.dimq();
  Eigen::VectorXd v = Eigen::VectorXd::Zero(fixed_base_robot_.dimv());
  Eigen::VectorXd vmax = fixed_base_robot_.jointVelocityLimit();
  ASSERT_EQ(dimq, fixed_base_robot_.jointVelocityLimit().size());
  limit.setSlackAndDual(dtau_, v);
  Eigen::VectorXd Cv = Eigen::VectorXd::Zero(dimq);
  limit.augmentDualResidual(dtau_, Cv);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimq);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq);
  slack_ref = dtau_ * (vmax-v);
  Eigen::VectorXd Cv_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::SetSlackAndDualPositive(dimq, barrier_, slack_ref, dual_ref);
  Cv_ref = dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  double cost_slack_barrier = limit.costSlackBarrier();
  double cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimq, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  Cv = Eigen::VectorXd::Ones(dimq);
  Eigen::MatrixXd Cvv = Eigen::MatrixXd::Ones(dimq, dimq);
  limit.condenseSlackAndDual(dtau_, v, Cvv, Cv);
  Cv_ref = Eigen::VectorXd::Ones(dimq);
  Eigen::MatrixXd Cvv_ref = Eigen::MatrixXd::Ones(dimq, dimq);
  for (int i=0; i<dimq; ++i) {
    Cvv_ref(i, i) += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  Eigen::VectorXd residual_ref = dtau_ * (v-vmax) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  Cv_ref.array() 
      += dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimq);
  limit.computeSlackAndDualDirection(dtau_, dv);
  const Eigen::VectorXd dslack_ref = - dtau_ * dv - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimq);
  pdipmfunc::ComputeDualDirection(dual_ref, slack_ref, dslack_ref, duality_ref, 
                                  ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(margin_rate);
  const double dual_step_size = limit.maxDualStepSize(margin_rate);
  const double slack_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimq, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimq, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(step_size);
  const double berrier_ref 
      = pdipmfunc::SlackBarrierCost(dimq, barrier_, 
                                    slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(step_size);
  limit.updateDual(step_size);
  cost_slack_barrier = limit.costSlackBarrier();
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimq, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  limit.augmentDualResidual(dtau_, Cv);
  Cv_ref += dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
}


TEST_F(JointVelocityUpperLimitsTest, condenseSlackAndDualFloatingBase) {
  JointVariablesUpperLimits limit(floating_base_robot_, 
                                  floating_base_robot_.jointVelocityLimit(), 
                                  barrier_);
  const int dimq = floating_base_robot_.dimq();
  const int dimv = floating_base_robot_.dimv();
  const int dim_passive = floating_base_robot_.dim_passive();
  Eigen::VectorXd v = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd vmax = floating_base_robot_.jointVelocityLimit();
  const int dimc = floating_base_robot_.jointVelocityLimit().size();
  ASSERT_EQ(dimc+6, dimv);
  limit.setSlackAndDual(dtau_, v);
  Eigen::VectorXd Cv = Eigen::VectorXd::Zero(dimv);
  limit.augmentDualResidual(dtau_, Cv);
  Eigen::VectorXd slack_ref = Eigen::VectorXd::Zero(dimc);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimc);
  slack_ref = dtau_ * (vmax-v.tail(dimc));
  Eigen::VectorXd Cv_ref = Eigen::VectorXd::Zero(dimv);
  pdipmfunc::SetSlackAndDualPositive(dimc, barrier_, slack_ref, dual_ref);
  Cv_ref.tail(dimc) = dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  double cost_slack_barrier = limit.costSlackBarrier();
  double cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimc, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  Cv = Eigen::VectorXd::Ones(dimv);
  Eigen::MatrixXd Cvv = Eigen::MatrixXd::Ones(dimv, dimv);
  limit.condenseSlackAndDual(dtau_, v, Cvv, Cv);
  Cv_ref = Eigen::VectorXd::Ones(dimv);
  Eigen::MatrixXd Cvv_ref = Eigen::MatrixXd::Ones(dimv, dimv);
  for (int i=0; i<dimc; ++i) {
    Cvv_ref(dim_passive+i, dim_passive+i) 
        += dtau_ * dtau_ * dual_ref.coeff(i) / slack_ref.coeff(i);
  }
  Eigen::VectorXd residual_ref = dtau_ * (v.tail(dimc)-vmax) + slack_ref;
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDualityResidual(barrier_, slack_ref, dual_ref, duality_ref);
  Cv_ref.tail(dimc).array() 
      += dtau_ * (dual_ref.array()*residual_ref.array()-duality_ref.array()) 
               / slack_ref.array();
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(dimv);
  limit.computeSlackAndDualDirection(dtau_, dv);
  const Eigen::VectorXd dslack_ref = - dtau_ * dv.tail(dimc) - residual_ref;
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dimc);
  pdipmfunc::ComputeDualDirection(dual_ref, slack_ref, dslack_ref, duality_ref, 
                                  ddual_ref);
  const double margin_rate = 0.995;
  const double slack_step_size = limit.maxSlackStepSize(margin_rate);
  const double dual_step_size = limit.maxDualStepSize(margin_rate);
  const double slack_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimc, margin_rate, slack_ref, dslack_ref);
  const double dual_step_size_ref 
      = pdipmfunc::FractionToBoundary(dimc, margin_rate, dual_ref, ddual_ref);
  EXPECT_DOUBLE_EQ(slack_step_size, slack_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
  const double step_size = std::min(slack_step_size, dual_step_size); 
  const double berrier = limit.costSlackBarrier(step_size);
  const double berrier_ref 
      = pdipmfunc::SlackBarrierCost(dimc, barrier_, 
                                    slack_ref+step_size*dslack_ref);
  EXPECT_DOUBLE_EQ(berrier, berrier_ref);
  limit.updateSlack(step_size);
  limit.updateDual(step_size);
  cost_slack_barrier = limit.costSlackBarrier();
  slack_ref += step_size * dslack_ref;
  dual_ref += step_size * ddual_ref;
  cost_slack_barrier_ref 
      = pdipmfunc::SlackBarrierCost(dimc, barrier_, slack_ref);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  limit.augmentDualResidual(dtau_, Cv);
  Cv_ref.tail(dimc) += dtau_ * dual_ref;
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
}

} // namespace pdipm
} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}