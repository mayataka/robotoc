#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {

class PDIPMTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dim_ = 100;
    barrier_ = 0.001;
    data_ = ConstraintComponentData(dim_);
    data_.slack.array() = Eigen::VectorXd::Random(dim_).array().abs();
    data_.dual.array() = Eigen::VectorXd::Random(dim_).array().abs();
    data_.dslack.array() = Eigen::VectorXd::Random(dim_);
    data_.ddual.array() = Eigen::VectorXd::Random(dim_);
  }

  virtual void TearDown() {
  }

  int dim_;
  double barrier_;
  ConstraintComponentData data_;
};


TEST_F(PDIPMTest, SetSlackAndDualPositive) {
  data_.slack = Eigen::VectorXd::Random(dim_);
  data_.dual = Eigen::VectorXd::Random(dim_);
  pdipm::SetSlackAndDualPositive(barrier_, data_);
  EXPECT_TRUE(data_.slack.minCoeff() >= barrier_);
  EXPECT_TRUE(data_.dual.minCoeff() >= barrier_);
}


TEST_F(PDIPMTest, ComputeDuality) {
  EXPECT_TRUE(data_.slack.minCoeff() >= 0);
  EXPECT_TRUE(data_.dual.minCoeff() >= 0);
  pdipm::ComputeDuality(barrier_, data_);
  Eigen::VectorXd duality_ref = Eigen::VectorXd::Zero(dim_);
  for (int i=0; i<dim_; ++i) {
    duality_ref(i) = data_.slack(i) * data_.dual(i) - barrier_;
  }
  EXPECT_TRUE(data_.duality.isApprox(duality_ref));
}


TEST_F(PDIPMTest, FractionToBoundary) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim_).array().abs();
  Eigen::VectorXd dvec = Eigen::VectorXd::Random(dim_);
  const double fraction_rate = 0.995;
  const double step_size = pdipm::FractionToBoundary(dim_, fraction_rate, 
                                                     vec, dvec);
  Eigen::VectorXd vec_updated = vec + step_size * dvec;
  EXPECT_TRUE(vec_updated.minCoeff() >= 0);
}


TEST_F(PDIPMTest, FractionToBoundarySlack) {
  EXPECT_TRUE(data_.slack.minCoeff() >= 0);
  const double fraction_rate = 0.995;
  const double step_slack = pdipm::FractionToBoundarySlack(fraction_rate, data_);
  Eigen::VectorXd slack_tmp = data_.slack + step_slack * data_.dslack;
  EXPECT_TRUE(slack_tmp.minCoeff() >= 0);
  const double step_size = pdipm::FractionToBoundary(dim_, fraction_rate, 
                                                     data_.slack, data_.dslack);
  EXPECT_DOUBLE_EQ(step_size, step_slack);
}


TEST_F(PDIPMTest, FractionToBoundaryDual) {
  EXPECT_TRUE(data_.dual.minCoeff() >= 0);
  const double fraction_rate = 0.995;
  const double step_dual = pdipm::FractionToBoundaryDual(fraction_rate, data_);
  Eigen::VectorXd dual_tmp = data_.dual + step_dual * data_.ddual;
  EXPECT_TRUE(dual_tmp.minCoeff() >= 0);
  const double step_size = pdipm::FractionToBoundary(dim_, fraction_rate, 
                                                     data_.dual, data_.ddual);
  EXPECT_DOUBLE_EQ(step_size, step_dual);
}



TEST_F(PDIPMTest, ComputeDualDirection) {
  data_.duality.array() = data_.dual.array() * data_.slack.array() - barrier_;
  pdipm::ComputeDualDirection(data_);
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dim_);
  for (int i=0; i<dim_; ++i) {
    ddual_ref(i) = - (data_.dual(i) * data_.dslack(i) + data_.duality(i)) / data_.slack(i);
  }
  EXPECT_TRUE(ddual_ref.isApprox(data_.ddual));
}


TEST_F(PDIPMTest, CostBarrier) {
  const double cost_ref = - barrier_ * (data_.slack.array().log()).sum();
  const double cost = pdipm::CostBarrier(barrier_, data_.slack);
  EXPECT_DOUBLE_EQ(cost_ref, cost);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}