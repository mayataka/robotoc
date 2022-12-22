#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/com_cost.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class CostFunctionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    t = grid_info.t;
    dt = grid_info.dt;
  }

  virtual void TearDown() {
  }

  void testStageCost(Robot& robot);

  GridInfo grid_info;
  double dt, t0, t;
};


void CostFunctionTest::testStageCost(Robot& robot) {
  const double discount_factor = 0.99;
  const double discount_time_step = 0.02;
  auto non_discounted_cost = std::make_shared<CostFunction>();
  auto discounted_cost = std::make_shared<CostFunction>(discount_factor, discount_time_step);
  EXPECT_DOUBLE_EQ(non_discounted_cost->discountFactor(), 1.0);
  EXPECT_DOUBLE_EQ(non_discounted_cost->discountTimeStep(), 0.0);
  EXPECT_DOUBLE_EQ(discounted_cost->discountFactor(), discount_factor);
  EXPECT_DOUBLE_EQ(discounted_cost->discountTimeStep(), discount_time_step);

  const int dimq = robot.dimq();
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  SplitKKTMatrix kkt_mat(robot);
  SplitKKTResidual kkt_res(robot);
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qaa.setRandom();
  kkt_mat.Quu.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu.setRandom();
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(dimv).array().abs().matrix();
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(dimv).array().abs().matrix(); 
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(dimv).array().abs().matrix();
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(dimu).array().abs().matrix();
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(dimv); 
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(dimu);
  auto config_cost = std::make_shared<ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(q_weight);
  config_cost->set_v_weight(v_weight);
  config_cost->set_a_weight(a_weight);
  config_cost->set_u_weight(u_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_ref(v_ref);
  config_cost->set_u_ref(u_ref);
  non_discounted_cost->add("config_cost", config_cost);
  discounted_cost->add("config_cost", config_cost);

  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();

  grid_info.stage = 10;

  auto data = CostFunctionData(robot);

  const auto s = SplitSolution::Random(robot, contact_status);

  const double non_discounted_value = non_discounted_cost->evalStageCost(robot, contact_status, grid_info, s, data);
  const double discounted_value = discounted_cost->evalStageCost(robot, contact_status, grid_info, s, data);
  // EXPECT_DOUBLE_EQ(non_discounted_value*std::pow(discount_factor, (grid_info.t-grid_info.t0)/discount_time_step), discounted_value);
  non_discounted_cost->setDiscountFactor(discount_factor, discount_time_step);
  EXPECT_DOUBLE_EQ(discounted_value, non_discounted_cost->evalStageCost(robot, contact_status, grid_info, s, data));

  auto cost_component = non_discounted_cost->get("config_cost");
  EXPECT_NO_THROW(
    auto config_cost_component = cost_component->as_shared_ptr<ConfigurationSpaceCost>();
    config_cost_component->set_q_weight(q_weight);
  );
  EXPECT_THROW(
    cost_component->as_shared_ptr<CoMCost>(),
    std::runtime_error
  );
  EXPECT_NO_THROW(
    std::cout << non_discounted_cost << std::endl;
    std::cout << *non_discounted_cost.get() << std::endl;
  );
}


TEST_F(CostFunctionTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  testStageCost(robot);
}


TEST_F(CostFunctionTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  testStageCost(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}