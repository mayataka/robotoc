#include <string>
#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/quadruped/cost_function.hpp"
#include "idocp/quadruped/constraints.hpp"


namespace idocp {

class FloatingBaseOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
    baum_on_velocity_ = 20;
    baum_on_position_ = 400;
    robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, baum_on_position_);
    for (int i=0; i<contact_frames_.size(); ++i) {
      contact_status_.push_back(rnd()%2==0);
    }
    robot_.setContactStatus(contact_status_);
    cost_ = std::make_shared<quadruped::CostFunction>(robot_);
    constraints_ = std::make_shared<quadruped::Constraints>(robot_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    N_ = 10;
    q_ = Eigen::VectorXd::Zero(robot_.dimq());
    q_ << 0, 0, 4.8, 0, 0, 0, 1, 
          0.0315, 0.4, -0.8, 
          0.0315, -0.4, 0.8, 
          -0.0315, 0.4, -0.8,
          -0.0315, -0.4, 0.8;
    v_ = Eigen::VectorXd::Zero(robot_.dimv());
    robot_.normalizeConfiguration(q_);
    robot_.updateKinematics(q_, v_, Eigen::VectorXd::Zero(robot_.dimv()));
    robot_.setContactPointsByCurrentKinematics();
    T_ = N_ * dtau_;
    ocp_ = OCP(robot_, cost_, constraints_, T_, N_);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  Robot robot_;
  std::shared_ptr<CostFunctionInterface> cost_;
  std::shared_ptr<ConstraintsInterface> constraints_;
  OCP ocp_;
  double t_, T_, dtau_, baum_on_velocity_, baum_on_position_;
  int N_, dimq_, dimv_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_;
  Eigen::VectorXd q_, v_;
  std::vector<int> contact_frames_;
  std::vector<bool> contact_status_;
  std::vector<std::vector<bool>> contact_sequence_;
};


TEST_F(FloatingBaseOCPTest, setStateTrajectory) {
  OCP ocp(robot_, cost_, constraints_, T_, N_);
  bool feasible = ocp.setStateTrajectory(q_, v_);
  EXPECT_TRUE(feasible);
  Eigen::VectorXd qN = Eigen::VectorXd::Zero(robot_.dimq());
  robot_.generateFeasibleConfiguration(qN);
  Eigen::VectorXd vN = Eigen::VectorXd::Zero(robot_.dimv());
  feasible = ocp.setStateTrajectory(q_, v_, qN, vN);
  EXPECT_TRUE(feasible);
}


TEST_F(FloatingBaseOCPTest, solveLQR) {
  bool use_line_search = false;
  ocp_.solveLQR(t_, q_, v_, use_line_search);
  use_line_search = true;
  ocp_.solveLQR(t_, q_, v_, use_line_search);
}


TEST_F(FloatingBaseOCPTest, KKTError) {
  bool use_line_search = false;
  std::cout << "contact status = [";
  for (int i=0; i<contact_status_.size(); ++i) {
    std::cout << contact_status_[i] << ", ";
  }
  std::cout << "]" << std::endl;;
  std::cout << "KKT error = " << ocp_.KKTError(t_, q_, v_) << std::endl;
  const int max_itr = 10;
  for (int i=0; i<max_itr; ++i) {
    ocp_.solveLQR(t_, q_, v_, use_line_search);
    std::cout << "KKT error = " << ocp_.KKTError(t_, q_, v_) << std::endl;
  }
  ocp_.printSolution();
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}