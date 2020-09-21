#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/ocp.hpp"


namespace idocp {

// tests the equivalence between the paralle and serial implementation.
class OCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    T_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    N_ = 10;
    num_proc_ = 1;
    dtau_ = T_ / N_;
  }

  virtual void TearDown() {
  }

  double dtau_, t_, T_;
  int N_, num_proc_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(OCPTest, updateSolutionFixedBaseWithoutContact) {
  Robot robot(fixed_base_urdf_);
  std::random_device rnd;
  auto cost = std::make_shared<CostFunction>();
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  auto contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.01);
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(robot.dimv());
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.01));
    f_ref.push_back(Eigen::Vector3d::Zero());
  }
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_vf_weight(vf_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
}


TEST_F(OCPTest, updateSolutionFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {rnd()%2==0};
  robot.setContactStatus(contact_status);
  auto cost = std::make_shared<CostFunction>();
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  auto contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Zero(robot.dimv());
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.01));
    f_ref.push_back(Eigen::Vector3d::Zero());
  }
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_vf_weight(vf_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  if (contact_status[0]) {
    ocp.activateContact(0, 0, N_);
  }
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  if (contact_status[0]) {
    ocp_ref.activateContact(0, 0, N_);
  }
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
}



TEST_F(OCPTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Constant(robot.dimv(), 10);
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Zero(robot.dimv());
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Constant(robot.dimv(), 0.1);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.01));
    f_ref.push_back(Eigen::Vector3d::Zero());
  }
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_vf_weight(vf_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  for (int i=0; i<contact_status.size(); ++i) {
    if (contact_status[i]) {
      ocp.activateContact(i, 0, N_);
    }
  }
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  for (int i=0; i<contact_status.size(); ++i) {
    if (contact_status[i]) {
      ocp_ref.activateContact(i, 0, N_);
    }
  }
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  EXPECT_DOUBLE_EQ(ocp.KKTError(t_), ocp_ref.KKTError(t_));
  EXPECT_DOUBLE_EQ(ocp.computeKKTError(t_, q, v), ocp_ref.computeKKTError(t_, q, v));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
