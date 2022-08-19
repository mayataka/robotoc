#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/constraints/pdipm.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"

namespace robotoc {

class PDIPMTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dim = 100;
    barrier_param = 0.001;
    data = ConstraintComponentData(dim, barrier_param);
    data.slack.array() = Eigen::VectorXd::Random(dim).array().abs();
    data.dual.array() = Eigen::VectorXd::Random(dim).array().abs();
    data.dslack.array() = Eigen::VectorXd::Random(dim);
    data.ddual.array() = Eigen::VectorXd::Random(dim);
  }

  virtual void TearDown() {
  }

  int dim;
  double barrier_param;
  ConstraintComponentData data;
};


TEST_F(PDIPMTest, setSlackAndDualPositive) {
  data.slack = Eigen::VectorXd::Random(dim);
  data.dual = Eigen::VectorXd::Random(dim);
  pdipm::setSlackAndDualPositive(barrier_param, data);
  EXPECT_TRUE(data.slack.minCoeff() >= barrier_param);
  EXPECT_TRUE(data.dual.minCoeff() >= barrier_param);
}


TEST_F(PDIPMTest, computeComplementarySlackness) {
  EXPECT_TRUE(data.slack.minCoeff() >= 0);
  EXPECT_TRUE(data.dual.minCoeff() >= 0);
  pdipm::computeComplementarySlackness(barrier_param, data);
  Eigen::VectorXd compl_ref = Eigen::VectorXd::Zero(dim);
  for (int i=0; i<dim; ++i) {
    compl_ref(i) = data.slack(i) * data.dual(i) - barrier_param;
  }
  EXPECT_TRUE(data.cmpl.isApprox(compl_ref));
}


TEST_F(PDIPMTest, fractionToBoundary) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(dim).array().abs();
  Eigen::VectorXd dvec = Eigen::VectorXd::Random(dim);
  const double fraction_rate = 0.995;
  const double step_size = pdipm::fractionToBoundary(dim, fraction_rate, 
                                                     vec, dvec);
  Eigen::VectorXd vec_updated = vec + step_size * dvec;
  EXPECT_TRUE(vec_updated.minCoeff() >= 0);
}


TEST_F(PDIPMTest, fractionToBoundarySlack) {
  EXPECT_TRUE(data.slack.minCoeff() >= 0);
  const double fraction_rate = 0.995;
  const double step_slack = pdipm::fractionToBoundarySlack(fraction_rate, data);
  Eigen::VectorXd slack_tmp = data.slack + step_slack * data.dslack;
  EXPECT_TRUE(slack_tmp.minCoeff() >= 0);
  const double step_size = pdipm::fractionToBoundary(dim, fraction_rate, 
                                                     data.slack, data.dslack);
  EXPECT_DOUBLE_EQ(step_size, step_slack);
}


TEST_F(PDIPMTest, fractionToBoundaryDual) {
  EXPECT_TRUE(data.dual.minCoeff() >= 0);
  const double fraction_rate = 0.995;
  const double step_dual = pdipm::fractionToBoundaryDual(fraction_rate, data);
  Eigen::VectorXd dual_tmp = data.dual + step_dual * data.ddual;
  EXPECT_TRUE(dual_tmp.minCoeff() >= 0);
  const double step_size = pdipm::fractionToBoundary(dim, fraction_rate, 
                                                     data.dual, data.ddual);
  EXPECT_DOUBLE_EQ(step_size, step_dual);
}



TEST_F(PDIPMTest, computeDualDirection) {
  data.cmpl.array() = data.dual.array() * data.slack.array() - barrier_param;
  pdipm::computeDualDirection(data);
  Eigen::VectorXd ddual_ref = Eigen::VectorXd::Zero(dim);
  for (int i=0; i<dim; ++i) {
    ddual_ref(i) = - (data.dual(i) * data.dslack(i) + data.cmpl(i)) / data.slack(i);
  }
  EXPECT_TRUE(ddual_ref.isApprox(data.ddual));
}


TEST_F(PDIPMTest, logBarrier) {
  const double cost_ref = - barrier_param * (data.slack.array().log()).sum();
  const double cost = pdipm::logBarrier(barrier_param, data.slack);
  EXPECT_DOUBLE_EQ(cost_ref, cost);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}