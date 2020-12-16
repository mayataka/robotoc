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
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/ocp.hpp"


namespace idocp {

class OCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);
  static std::shared_ptr<Constraints> createConstraints(const Robot& robot);
  static void testUpdateSolutionInParallel(Robot& robot);
  static void testUpdateSolutionInParallelWithoutActiveContacts(Robot& robot);
  static void testUpdateSolutionInParallelWithActiveContacts(Robot& robot);
  std::string fixed_base_urdf, floating_base_urdf;
};


std::shared_ptr<CostFunction> OCPTest::createCost(const Robot& robot) {
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimu()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_a_ref(a_ref);
  joint_cost->set_u_weight(u_weight);
  joint_cost->set_u_ref(u_ref);
  joint_cost->set_qf_weight(qf_weight);
  joint_cost->set_vf_weight(vf_weight);
  const int task_frame = 10;
  auto contact_force_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_force_cost->set_f_weight(f_weight);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(joint_cost);
  cost->push_back(contact_force_cost);
  return cost;
}


std::shared_ptr<Constraints> OCPTest::createConstraints(const Robot& robot) {
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto torques_lower_limit = std::make_shared<JointTorquesLowerLimit>(robot);
  auto torques_upper_limit = std::make_shared<JointTorquesUpperLimit>(robot);
  auto constraints = std::make_shared<Constraints>();
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  constraints->push_back(torques_lower_limit);
  constraints->push_back(torques_upper_limit);
  return constraints;
}


void OCPTest::testUpdateSolutionInParallel(Robot& robot) {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
  const double t = 0;
  const double T = 1;
  const double N = 8;
  Eigen::VectorXd q(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T, N, 1);
  OCP ocp_ref(robot, cost, constraints, T, N, 4);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, false);
  ocp_ref.updateSolution(t, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, true);
  ocp_ref.updateSolution(t, q, v, true);
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(ocp.getSolution(i).isApprox(ocp_ref.getSolution(i)));
  }
}


void OCPTest::testUpdateSolutionInParallelWithoutActiveContacts(Robot& robot) {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
  const double t = 0;
  const double T = 1;
  const double N = 8;
  Eigen::VectorXd q(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T, N, 1);
  OCP ocp_ref(robot, cost, constraints, T, N, 4);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    ocp.deactivateContacts({i}, 0, N);
    ocp_ref.deactivateContacts({i}, 0, N);
  }
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, false);
  ocp_ref.updateSolution(t, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, true);
  ocp_ref.updateSolution(t, q, v, true);
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(ocp.getSolution(i).isApprox(ocp_ref.getSolution(i)));
  }
}


void OCPTest::testUpdateSolutionInParallelWithActiveContacts(Robot& robot) {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
  const double t = 0;
  const double T = 1;
  const double N = 8;
  Eigen::VectorXd q(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  OCP ocp(robot, cost, constraints, T, N, 1);
  OCP ocp_ref(robot, cost, constraints, T, N, 4);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    ocp.activateContacts({i}, 0, N);
    ocp_ref.activateContacts({i}, 0, N);
  }
  ocp.setStateTrajectory(q, v);
  ocp_ref.setStateTrajectory(q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, false);
  ocp_ref.updateSolution(t, q, v, false);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  ocp.updateSolution(t, q, v, true);
  ocp_ref.updateSolution(t, q, v, true);
  ocp.computeKKTResidual(t, q, v);
  ocp_ref.computeKKTResidual(t, q, v);
  EXPECT_DOUBLE_EQ(ocp.KKTError(), ocp_ref.KKTError());
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(ocp.getSolution(i).isApprox(ocp_ref.getSolution(i)));
  }
}


TEST_F(OCPTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testUpdateSolutionInParallel(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  testUpdateSolutionInParallelWithoutActiveContacts(robot);
  testUpdateSolutionInParallelWithActiveContacts(robot);
}


TEST_F(OCPTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testUpdateSolutionInParallel(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  testUpdateSolutionInParallelWithoutActiveContacts(robot);
  testUpdateSolutionInParallelWithActiveContacts(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
