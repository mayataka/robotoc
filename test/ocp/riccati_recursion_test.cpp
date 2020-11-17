#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_recursion.hpp"

namespace idocp {

class RiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    N = 20;
    max_num_impulse = 5;
    nproc = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  std::shared_ptr<CostFunctionData> createCost() const;
  std::shared_ptr<Constraints> createConstraints() const;
  void testBackwardRiccatiRecursionTerminal(const Robot& robot_const) const;
  void testBackwardRiccatiRecursion(const Robot& robot_const) const;

  int N, max_num_impulse, nproc;
  double T, t, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


std::shared_ptr<CostFunctionData> RiccatiRecursionTest::createCost(const Robot& robot) const {
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
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_force_cost->set_f_weight(f_weight);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(joint_cost);
  cost->push_back(contact_force_cost);
  return cost;
}


std::shared_ptr<Constraints> RiccatiRecursionTest::createConstraints(const Robot& robot) const {
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


void RiccatiRecursionTest::testBackwardRiccatiRecursionTerminal(const Robot& robot_const) const {
  TerminalOCP terminal_ocp(robot createCost(robot_const), createConstraints(robot_const));
  auto robot = robot_const;
  const auto s = SplitSolution::Random(robot_const);
  terminal_ocp.linearizeOCP(robot, t, s);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(terminal_ocp);

  auto robot_ref = robot_const;
  TerminalOCP terminal_ocp_ref(robot createCost(robot_const), createConstraints(robot_const));
  terminal_ocp_ref.linearizeOCP(robot_ref, t, s);
  
}


TEST_F(RiccatiRecursionTest, fixedBase) {
  testBackwardRecursion(fixed_base_robot);
  testForwardRecursionWithoutStateConstraint(fixed_base_robot);
  testForwardRecursionWithStateConstraint(fixed_base_robot);
  testFactorizeStateConstraintFactorization(fixed_base_robot);
}


TEST_F(RiccatiRecursionTest, floating_base) {
  testBackwardRecursion(floating_base_robot);
  testForwardRecursionWithoutStateConstraint(floating_base_robot);
  testForwardRecursionWithStateConstraint(floating_base_robot);
  testFactorizeStateConstraintFactorization(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}