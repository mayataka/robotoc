#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/ocp.hpp"


namespace idocp {

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
  // auto contact_cost = std::make_shared<ContactCost>(robot);
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
  // contact_cost->set_f_weight(f_weight);
  // contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  // cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  constraints->push_back(torques_lower_limit); 
  constraints->push_back(torques_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::cout << "Aaa" << std::endl;
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  std::cout << "Aaa" << std::endl;
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
}


TEST_F(OCPTest, updateSolutionFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  auto cost = std::make_shared<CostFunction>();
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  // auto contact_cost = std::make_shared<ContactCost>(robot);
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
  // contact_cost->set_f_weight(f_weight);
  // contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  // cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  constraints->push_back(torques_lower_limit); 
  constraints->push_back(torques_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  std::random_device rnd;
  // const bool is_contact_active = (rnd()%2==0);
  const bool is_contact_active = true;
  if (is_contact_active) {
    ocp.activateContact(0, 0, N_);
  }
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  if (is_contact_active) {
    ocp_ref.activateContact(0, 0, N_);
  }
  ocp.setStateTrajectory(q, v);
  ocp_ref.setStateTrajectory(q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
}


TEST_F(OCPTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::shared_ptr<CostFunction> cost = std::make_shared<CostFunction>();
  std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
  // std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
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
  // contact_cost->set_f_weight(f_weight);
  // contact_cost->set_f_ref(f_ref);
  cost->push_back(joint_cost);
  // cost->push_back(contact_cost);
  auto constraints = std::make_shared<Constraints>();
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  constraints->push_back(torques_lower_limit); 
  constraints->push_back(torques_upper_limit);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T_, N_, 1);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (int i=0; i<contact_frames.size(); ++i) {
    is_contact_active.push_back(rnd()%2==0);
  }
  for (int i=0; i<contact_frames.size(); ++i) {
    if (is_contact_active[i]) {
      ocp.activateContact(i, 0, N_);
    }
  }
  OCP ocp_ref(robot, cost, constraints, T_, N_, 2);
  for (int i=0; i<contact_frames.size(); ++i) {
    if (is_contact_active[i]) {
      ocp_ref.activateContact(i, 0, N_);
    }
  }
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, false);
  ocp_ref.updateSolution(t_, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t_, q, v, true);
  ocp_ref.updateSolution(t_, q, v, true);
  ocp.computeKKTResidual(t_, q, v);
  ocp_ref.computeKKTResidual(t_, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
