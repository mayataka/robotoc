#include <string>
#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/manipulator/cost_function.hpp"
#include "idocp/manipulator/constraints.hpp"


namespace idocp {

class FixedBaseOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    contact_frames_ = {18};
    baum_on_velocity_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    baum_on_position_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    if (rnd()%2==0) {
      robot_ = Robot(urdf_, contact_frames_, baum_on_velocity_, 
                     baum_on_position_);
    }
    else {
      robot_ = Robot(urdf_);
    }
    cost_ = std::make_shared<manipulator::CostFunction>(robot_);
    constraints_ = std::make_shared<manipulator::Constraints>(robot_);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    T_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    N_ = 10;
    dtau_ = T_/N_;
    if (rnd()%2==0) {
      contact_sequence_ = std::vector<std::vector<bool>>(N_, {true});
    }
    else {
      contact_sequence_ = std::vector<std::vector<bool>>(N_, {false});
    }
    q_ = Eigen::VectorXd::Zero(robot_.dimq());
    robot_.generateFeasibleConfiguration(q_);
    v_ = Eigen::VectorXd::Random(robot_.dimv());
    if (robot_.max_point_contacts() > 0) {
      robot_.updateKinematics(q_, v_, Eigen::VectorXd::Zero(robot_.dimv()));
      robot_.setContactPointsByCurrentKinematics();
    }
    std::cout << "T = " << T_ << std::endl;
    ocp_ = OCP(robot_, cost_, constraints_, T_, N_);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  Robot robot_;
  std::shared_ptr<CostFunctionInterface> cost_;
  std::shared_ptr<ConstraintsInterface> constraints_;
  OCP ocp_;
  std::vector<int> contact_frames_;
  std::vector<bool> contact_status_;
  std::vector<std::vector<bool>> contact_sequence_;
  double t_, T_, dtau_, baum_on_velocity_, baum_on_position_;
  int N_, dimq_, dimv_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_;
  Eigen::VectorXd q_, v_;
};


TEST_F(FixedBaseOCPTest, setStateTrajectory) {
  OCP ocp(robot_, cost_, constraints_, T_, N_);
  bool feasible = ocp.setStateTrajectory(q_, v_);
  EXPECT_TRUE(feasible);
  Eigen::VectorXd qN = Eigen::VectorXd::Zero(robot_.dimq());
  robot_.generateFeasibleConfiguration(qN);
  Eigen::VectorXd vN = Eigen::VectorXd::Zero(robot_.dimv());
  feasible = ocp.setStateTrajectory(q_, v_, qN, vN);
  EXPECT_TRUE(feasible);
}


TEST_F(FixedBaseOCPTest, solveLQR) {
  std::cout << "contact is ";
  if (robot_.max_dimf() == 0) {
    std::cout << "not set in the model" << std::endl;
  }
  else if (contact_sequence_[0][0]) {
    ocp_.setContactSequence(contact_sequence_);
    std::cout << "active" << std::endl;
  }
  else {
    ocp_.setContactSequence(contact_sequence_);
    std::cout << "not active" << std::endl;
  }
  bool use_line_search = false;
  ocp_.solveLQR(t_, q_, v_, use_line_search);
  use_line_search = true;
  ocp_.solveLQR(t_, q_, v_, use_line_search);
}


TEST_F(FixedBaseOCPTest, KKTError) {
  std::cout << "contact is ";
  if (robot_.max_dimf() == 0) {
    std::cout << "not set in the model" << std::endl;
  }
  else if (contact_sequence_[0][0]) {
    ocp_.setContactSequence(contact_sequence_);
    std::cout << "active" << std::endl;
  }
  else {
    ocp_.setContactSequence(contact_sequence_);
    std::cout << "not active" << std::endl;
  }
  std::cout << "KKT error = " << ocp_.KKTError(t_, q_, v_) << std::endl;
  ocp_.solveLQR(t_, q_, v_);
  std::cout << "KKT error = " << ocp_.KKTError(t_, q_, v_) << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}