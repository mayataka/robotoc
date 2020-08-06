#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/quadruped/cost_function.hpp"


namespace idocp {

class QuadrupedCostFunctionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
    robot_ = Robot(urdf_, contact_frames_, 0, 0);
    std::vector<bool> contact_status;
    for (int i=0; i<contact_frames_.size(); ++i) {
      contact_status.push_back(rnd()%2==0);
    }
    robot_.setContactStatus(contact_status);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_ = Eigen::VectorXd::Zero(robot_.dimq());
    robot_.generateFeasibleConfiguration(q_);
    v_ = Eigen::VectorXd::Random(robot_.dimv());
    a_ = Eigen::VectorXd::Random(robot_.dimv());
    u_ = Eigen::VectorXd::Random(robot_.dimv());
    f_ = Eigen::VectorXd::Random(robot_.max_dimf());
    data_ = CostFunctionData(robot_);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
  CostFunctionData data_;
  Eigen::VectorXd q_, v_, a_, u_, f_;
};


TEST_F(QuadrupedCostFunctionTest, moveAssign) {
  std::unique_ptr<CostFunctionInterface> cost 
      = std::make_unique<quadruped::CostFunction>(robot_);
  std::unique_ptr<CostFunctionInterface> cost_ref;
  const double l = cost->l(robot_, data_, 0, dtau_, q_, v_, a_, u_, f_);
  cost_ref = std::move(cost);
  const double l_ref = cost_ref->l(robot_, data_, 0, dtau_, q_, v_, a_, u_, f_);
  EXPECT_DOUBLE_EQ(l, l_ref);
}


TEST_F(QuadrupedCostFunctionTest, moveConstructor) {
  std::unique_ptr<CostFunctionInterface> cost 
      = std::make_unique<quadruped::CostFunction>(robot_);
  const double l = cost->l(robot_, data_, 0, dtau_, q_, v_, a_, u_, f_);
  std::unique_ptr<CostFunctionInterface> cost_ref(std::move(cost));
  const double l_ref = cost_ref->l(robot_, data_, 0, dtau_, q_, v_, a_, u_, f_);
  EXPECT_DOUBLE_EQ(l, l_ref);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}