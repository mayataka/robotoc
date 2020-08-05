#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/manipulator/cost_function.hpp"


namespace idocp {

class ManipulatorCostFunctionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    contact_frames_ = {18};
    robot_ = Robot(urdf_, contact_frames_, 0, 0);
    std::vector<bool> contact_status = {rnd()%2==0};
    robot_.setContactStatus(contact_status);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_ = Eigen::VectorXd::Zero(robot_.dimq());
    robot_.generateFeasibleConfiguration(q_);
    v_ = Eigen::VectorXd::Random(robot_.dimv());
    a_ = Eigen::VectorXd::Random(robot_.dimv());
    u_ = Eigen::VectorXd::Random(robot_.dimv());
    f_ = Eigen::VectorXd::Random(robot_.max_dimf());
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
  Eigen::VectorXd q_, v_, a_, u_, f_;
};


TEST_F(ManipulatorCostFunctionTest, moveAssign) {
  std::unique_ptr<CostFunctionInterface> cost 
      = std::make_unique<manipulator::CostFunction>(robot_);
  std::unique_ptr<CostFunctionInterface> cost_ref;
  const double l = cost->l(robot_, 0, dtau_, q_, v_, a_, u_, f_);
  cost_ref = std::move(cost);
  const double l_ref = cost_ref->l(robot_, 0, dtau_, q_, v_, a_, u_, f_);
  EXPECT_DOUBLE_EQ(l, l_ref);
}


TEST_F(ManipulatorCostFunctionTest, moveConstructor) {
  std::unique_ptr<CostFunctionInterface> cost 
      = std::make_unique<manipulator::CostFunction>(robot_);
  const double l = cost->l(robot_, 0, dtau_, q_, v_, a_, u_, f_);
  std::unique_ptr<CostFunctionInterface> cost_ref(std::move(cost));
  const double l_ref = cost_ref->l(robot_, 0, dtau_, q_, v_, a_, u_, f_);
  EXPECT_DOUBLE_EQ(l, l_ref);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}