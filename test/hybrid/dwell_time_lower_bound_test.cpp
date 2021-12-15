#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/hybrid/time_discretization.hpp"
#include "robotoc/hybrid/dwell_time_lower_bound.hpp"


namespace robotoc {

class DwellTimeLowerBoundTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    min_dt = Eigen::VectorXd::Random(1).array().abs()[0];
    barrier = 1e-03;
    fraction_to_boundary_rule = 0.995;
    bound = DwellTimeLowerBound(barrier, fraction_to_boundary_rule);
    ts1 = Eigen::VectorXd::Random(1).array().abs()[0];
    ts2 = ts1 + Eigen::VectorXd::Random(1).array().abs()[0];
  }

  virtual void TearDown() {
  }

  int max_point_contacts;

  DwellTimeLowerBound bound;
  double min_dt, barrier, fraction_to_boundary_rule, ts1, ts2;
};


TEST_F(DwellTimeLowerBoundTest, testParam) {
  EXPECT_DOUBLE_EQ(bound.barrier(), 1.0e-03);
  EXPECT_DOUBLE_EQ(bound.fractionToBoundaryRule(), 0.995);
  // overwrite parameters
  const double barrier = 0.8;
  const double fraction_to_boundary_rule = 0.5;
  bound.setBarrier(barrier);
  bound.setFractionToBoundaryRule(fraction_to_boundary_rule);
  EXPECT_DOUBLE_EQ(bound.barrier(), barrier);
  EXPECT_DOUBLE_EQ(bound.fractionToBoundaryRule(), fraction_to_boundary_rule);
}


TEST_F(DwellTimeLowerBoundTest, lub) {
  // overwrite parameters
  const double barrier = 1e-03;
  const double fraction_to_boundary_rule = 0.95;
  bound.setBarrier(barrier);
  bound.setFractionToBoundaryRule(fraction_to_boundary_rule);
  bound.setSlack(min_dt, ts1, ts2);
  double slack_ref = ts2 - ts1 - min_dt;
  if (slack_ref <= 0) {
    slack_ref = std::sqrt(barrier);
  }
  double dual_ref = barrier / slack_ref;

  bound.evalConstraint(min_dt, ts1, ts2);
  double residual_ref = ts1 + min_dt - ts2 + slack_ref;
  double cmpl_ref = slack_ref * dual_ref - barrier;
  const double log_barrier = - barrier * std::log(slack_ref);
  const double kkt_error = bound.KKTError();
  const double kkt_error_ref = residual_ref*residual_ref + cmpl_ref*cmpl_ref;
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  SplitKKTResidual kkt_res1, kkt_res2;
  bound.evalDerivatives_lub(kkt_res1, kkt_res2);
  double h1_ref = dual_ref;
  double h2_ref = - dual_ref;
  EXPECT_DOUBLE_EQ(kkt_res1.h, -h1_ref);
  EXPECT_DOUBLE_EQ(kkt_res2.h, -h2_ref);
  SplitKKTMatrix kkt_mat1, kkt_mat2;
  bound.condenseSlackAndDual_lub(kkt_mat1, kkt_res1, kkt_mat2, kkt_res2);
  double Qtt1_ref = dual_ref / slack_ref;
  h1_ref += (dual_ref*residual_ref-cmpl_ref) / slack_ref;
  double Qtt2_ref = dual_ref / slack_ref;
  h2_ref -= (dual_ref*residual_ref-cmpl_ref) / slack_ref;
  EXPECT_DOUBLE_EQ(kkt_mat1.Qtt, Qtt1_ref);
  EXPECT_DOUBLE_EQ(kkt_res1.h, -h1_ref);
  EXPECT_DOUBLE_EQ(kkt_mat2.Qtt, Qtt2_ref);
  EXPECT_DOUBLE_EQ(kkt_res2.h, -h2_ref);
}


TEST_F(DwellTimeLowerBoundTest, lb) {
  // In this case, ts1 is a fixed parameter and ts2 is an optimization variable.
  bound.setSlack(min_dt, ts1, ts2);
  double slack_ref = ts2 - ts1 - min_dt;
  if (slack_ref <= 0) {
    slack_ref = std::sqrt(barrier);
  }
  double dual_ref = barrier / slack_ref;

  bound.evalConstraint(min_dt, ts1, ts2);
  double residual_ref = ts1 + min_dt - ts2 + slack_ref;
  double cmpl_ref = slack_ref * dual_ref - barrier;
  const double log_barrier = - barrier * std::log(slack_ref);
  const double kkt_error = bound.KKTError();
  const double kkt_error_ref = residual_ref*residual_ref + cmpl_ref*cmpl_ref;
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  SplitKKTResidual kkt_res2;
  bound.evalDerivatives_lb(kkt_res2);
  double h2_ref = - dual_ref;
  EXPECT_DOUBLE_EQ(kkt_res2.h, -h2_ref);
  SplitKKTMatrix kkt_mat2;
  bound.condenseSlackAndDual_lb(kkt_mat2, kkt_res2);
  double Qtt2_ref = dual_ref / slack_ref;
  h2_ref -= (dual_ref*residual_ref-cmpl_ref) / slack_ref;
  EXPECT_DOUBLE_EQ(kkt_mat2.Qtt, Qtt2_ref);
  EXPECT_DOUBLE_EQ(kkt_res2.h, -h2_ref);
}


TEST_F(DwellTimeLowerBoundTest, ub) {
  // In this case, ts2 is a fixed parameter and ts1 is an optimization variable.
  bound.setSlack(min_dt, ts1, ts2);
  double slack_ref = ts2 - ts1 - min_dt;
  if (slack_ref <= 0) {
    slack_ref = std::sqrt(barrier);
  }
  double dual_ref = barrier / slack_ref;

  bound.evalConstraint(min_dt, ts1, ts2);
  double residual_ref = ts1 + min_dt - ts2 + slack_ref;
  double cmpl_ref = slack_ref * dual_ref - barrier;
  const double log_barrier = - barrier * std::log(slack_ref);
  const double kkt_error = bound.KKTError();
  const double kkt_error_ref = residual_ref*residual_ref + cmpl_ref*cmpl_ref;
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  SplitKKTResidual kkt_res1;
  bound.evalDerivatives_ub(kkt_res1);
  double h1_ref = dual_ref;
  EXPECT_DOUBLE_EQ(kkt_res1.h, -h1_ref);
  SplitKKTMatrix kkt_mat1;
  bound.condenseSlackAndDual_ub(kkt_mat1, kkt_res1);
  double Qtt1_ref = dual_ref / slack_ref;
  h1_ref += (dual_ref*residual_ref-cmpl_ref) / slack_ref;
  EXPECT_DOUBLE_EQ(kkt_mat1.Qtt, Qtt1_ref);
  EXPECT_DOUBLE_EQ(kkt_res1.h, -h1_ref);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}