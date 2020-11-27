#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/ocp_riccati_solver.hpp"


namespace idocp {

class OCPRiccatiSolverTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_impulse = 5;
    nproc = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);
  static std::shared_ptr<Constraints> createConstraints(const Robot& robot);
  HybridSolution createSolution(const Robot& robot) const;
  HybridSolution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  template <typename T>
  void testIsSame(const T& rhs, const T& lhs) const {
    for (int i=0; i<=N; ++i) {
      EXPECT_TRUE(rhs[i].isApprox(lhs[i]));
      EXPECT_FALSE(rhs[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.impulse[i].isApprox(lhs.impulse[i]));
      EXPECT_FALSE(rhs.impulse[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.aux[i].isApprox(lhs.aux[i]));
      EXPECT_FALSE(rhs.aux[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.lift[i].isApprox(lhs.lift[i]));
      EXPECT_FALSE(rhs.lift[i].hasNaN());
    }
  }

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, t, dtau;
};


std::shared_ptr<CostFunction> OCPRiccatiSolverTest::createCost(const Robot& robot) {
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
  auto impulse_joint_cost = std::make_shared<JointSpaceImpulseCost>(robot);
  impulse_joint_cost->set_q_weight(q_weight);
  impulse_joint_cost->set_q_ref(q_ref);
  impulse_joint_cost->set_v_weight(v_weight);
  impulse_joint_cost->set_v_ref(v_ref);
  impulse_joint_cost->set_dv_weight(a_weight);
  impulse_joint_cost->set_dv_ref(a_ref);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  cost->push_back(joint_cost);
  cost->push_back(contact_force_cost);
  cost->push_back(impulse_joint_cost);
  cost->push_back(impulse_force_cost);
  return cost;
}


std::shared_ptr<Constraints> OCPRiccatiSolverTest::createConstraints(const Robot& robot) {
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


HybridSolution OCPRiccatiSolverTest::createSolution(const Robot& robot) const {
  HybridSolution s(N+1, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


HybridSolution OCPRiccatiSolverTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  if (robot.max_point_contacts() == 0) {
    return createSolution(robot);
  }
  else {
    HybridSolution s(N+1, max_num_impulse, robot);
    for (int i=0; i<=N; ++i) {
      s[i].setRandom(robot, contact_sequence.contactStatus(i));
    }
    const int num_impulse = contact_sequence.totalNumImpulseStages();
    for (int i=0; i<num_impulse; ++i) {
      s.impulse[i].setRandom(robot, contact_sequence.impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      const int time_stage = contact_sequence.timeStageBeforeImpulse(i);
      s.aux[i].setRandom(robot, contact_sequence.contactStatus(time_stage+1));
    }
    const int num_lift = contact_sequence.totalNumLiftStages();
    for (int i=0; i<num_lift; ++i) {
      const int time_stage = contact_sequence.timeStageBeforeLift(i);
      s.lift[i].setRandom(robot, contact_sequence.contactStatus(time_stage+1));
    }
    return s;
  }
}


ContactSequence OCPRiccatiSolverTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, T, N);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  for (int i=0; i<max_num_impulse; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    tmp.eventTime = i * 0.15 + 0.1 * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void OCPRiccatiSolverTest::test(const Robot& robot) const {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
  OCPLinearizer linearizer(T, N, max_num_impulse, nproc);
  ContactSequence contact_sequence(robot, T, N);
  if (robot.max_point_contacts() > 0) {
    contact_sequence = createContactSequence(robot);
  }
  auto kkt_matrix = HybridKKTMatrix(N+1, max_num_impulse, robot);
  auto kkt_residual = HybridKKTResidual(N+1, max_num_impulse, robot);
  HybridSolution s;
  if (robot.max_point_contacts() > 0) {
    s = createSolution(robot, contact_sequence);
  }
  else {
    s = createSolution(robot);
  }
  Eigen::VectorXd q(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto split_ocps = HybridOCP(N, max_num_impulse, robot, cost, constraints);
  std::vector<Robot> robots(nproc, robot);
  linearizer.initConstraints(split_ocps, robots, contact_sequence, s);
  linearizer.linearizeOCP(split_ocps, robots, contact_sequence, t, q, v, s, kkt_matrix, kkt_residual);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  auto split_ocps_ref = split_ocps;
  HybridDirection d = HybridDirection(N, max_num_impulse, robot);
  for (int i=0; i<N; ++i) {
    d[i].setContactStatus(contact_sequence.contactStatus(i));
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<num_impulse; ++i) {
    d.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d.aux[i].setContactStatus(contact_sequence.contactStatus(contact_sequence.timeStageBeforeImpulse(i)+1));
  }
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_lift; ++i) {
    d.lift[i].setContactStatus(contact_sequence.contactStatus(contact_sequence.timeStageBeforeLift(i)+1));
  }
  auto d_ref = d;
  OCPRiccatiSolver riccati_solver(robot, T, N, max_num_impulse, nproc);
  riccati_solver.setConstraintDimensions(contact_sequence);
  riccati_solver.computeDirection(split_ocps, robots, contact_sequence, q, v, s, d, kkt_matrix, kkt_residual);

}


TEST_F(OCPRiccatiSolverTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  test(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(OCPRiccatiSolverTest, floatingBase) {
  Robot robot(floating_base_urdf);
  test(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}