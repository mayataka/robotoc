#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/quadruped/constraints.hpp"


namespace idocp {

class QuadrupedConstraintsTest : public ::testing::Test {
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
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
  Eigen::VectorXd q_, v_, a_, u_, f_;
};


TEST_F(QuadrupedConstraintsTest, moveAssign) {
  std::unique_ptr<ConstraintsInterface> constraints
      = std::make_unique<quadruped::Constraints>(robot_);
  std::unique_ptr<ConstraintsInterface> constraints_ref;
  constraints_ref = std::move(constraints);
}


TEST_F(QuadrupedConstraintsTest, moveConstructor) {
  std::unique_ptr<ConstraintsInterface> constraints
      = std::make_unique<quadruped::Constraints>(robot_);
  std::unique_ptr<ConstraintsInterface> constraints_ref(std::move(constraints));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}