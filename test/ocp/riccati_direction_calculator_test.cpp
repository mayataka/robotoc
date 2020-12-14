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
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/riccati_direction_calculator.hpp"


namespace idocp {

class RiccatiDirectionCalculatorTest : public ::testing::Test {
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
  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  template <typename T>
  void testIsSame(const T& rhs, const T& lhs) const {
    for (int i=0; i<=N; ++i) {
      EXPECT_TRUE(rhs[i].isApprox(lhs[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.impulse[i].isApprox(lhs.impulse[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.aux[i].isApprox(lhs.aux[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.lift[i].isApprox(lhs.lift[i]));
    }
  }

  template <typename T>
  void testNaN(const T& ob) const {
    for (int i=0; i<=N; ++i) {
      EXPECT_FALSE(ob[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_FALSE(ob.impulse[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_FALSE(ob.aux[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_FALSE(ob.lift[i].hasNaN());
    }
  }

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, t, dtau;
};


std::shared_ptr<CostFunction> RiccatiDirectionCalculatorTest::createCost(const Robot& robot) {
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


std::shared_ptr<Constraints> RiccatiDirectionCalculatorTest::createConstraints(const Robot& robot) {
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


Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot) const {
  Solution s(N, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  if (robot.max_point_contacts() == 0) {
    return createSolution(robot);
  }
  else {
    Solution s(N, max_num_impulse, robot);
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


ContactSequence RiccatiDirectionCalculatorTest::createContactSequence(const Robot& robot) const {
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
    tmp.eventTime = i * 0.15 + 0.01 * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void RiccatiDirectionCalculatorTest::test(const Robot& robot) const {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
  ContactSequence contact_sequence(robot, T, N);
  if (robot.max_point_contacts() > 0) {
    contact_sequence = createContactSequence(robot);
  }
  auto kkt_matrix = KKTMatrix(N, max_num_impulse, robot);
  auto kkt_residual = KKTResidual(N, max_num_impulse, robot);
  Solution s;
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
  OCPLinearizer linearizer(T, N, max_num_impulse, nproc);
  std::vector<Robot> robots(nproc, robot);
  linearizer.initConstraints(split_ocps, robots, contact_sequence, s);
  linearizer.linearizeOCP(split_ocps, robots, contact_sequence, t, q, v, s, kkt_matrix, kkt_residual);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  StateConstraintRiccatiFactorizer constraint_factorizer(robot, N, max_num_impulse, nproc);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, 
                                             kkt_matrix, kkt_residual, 
                                             factorization);
  if (contact_sequence.existImpulseStage()) {
    riccati_recursion.forwardStateConstraintFactorization(
        factorizer, contact_sequence, kkt_matrix, kkt_residual, 
        factorization);
    riccati_recursion.backwardStateConstraintFactorization(
        factorizer, contact_sequence, kkt_matrix, 
        constraint_factorization);
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  const int num_lift = contact_sequence.totalNumLiftStages();
  Direction d = Direction(N, max_num_impulse, robot);
  for (int i=0; i<N; ++i) {
    d[i].setContactStatus(contact_sequence.contactStatus(i));
  }
  for (int i=0; i<num_impulse; ++i) {
    d.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d.impulse[i].dxi().setRandom();
    d.aux[i].setContactStatus(contact_sequence.contactStatus(contact_sequence.timeStageBeforeImpulse(i)+1));
  }
  for (int i=0; i<num_lift; ++i) {
    d.lift[i].setContactStatus(contact_sequence.contactStatus(contact_sequence.timeStageBeforeLift(i)+1));
  }
  auto d_ref = d;
  RiccatiDirectionCalculator::computeInitialStateDirection(robots, q, v, s, d);
  robot.subtractConfiguration(q, s[0].q, d_ref[0].dq());
  d_ref[0].dv() = v - s[0].v;
  testIsSame(d_ref, d);
  if (contact_sequence.existImpulseStage()) {
    constraint_factorizer.computeLagrangeMultiplierDirection(
        contact_sequence, factorization, constraint_factorization, d);
    constraint_factorizer.aggregateLagrangeMultiplierDirection(
        constraint_factorization, contact_sequence, d, factorization);
  }
  riccati_recursion.forwardRiccatiRecursion(
      factorizer, contact_sequence, kkt_matrix, kkt_residual, 
      factorization, d);
  testNaN(factorization);
  testNaN(kkt_matrix);
  testNaN(kkt_residual);
  d_ref = d;
  auto split_ocps_ref = split_ocps;
  RiccatiDirectionCalculator direction_calculator(T, N, max_num_impulse, nproc);
  direction_calculator.computeNewtonDirectionFromRiccatiFactorization(
      split_ocps, robots, contact_sequence, factorizer, 
      factorization, s, d);
  const double primal_step_size = direction_calculator.maxPrimalStepSize();
  const double dual_step_size = direction_calculator.maxDualStepSize();
  const Eigen::VectorXd dx0 = d_ref[0].dx();
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  double primal_step_size_ref = 1;
  double dual_step_size_ref = 1;
  auto robot_ref = robots[0];
  for (int i=0; i<N; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double impulse_time = contact_sequence.impulseTime(impulse_index);
      const double dtau_impulse = impulse_time - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                 exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization.impulse[impulse_index], d_ref[i], 
          exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(robot_ref, dtau_impulse,
                                                        s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
      // impulse
      ImpulseSplitRiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d_ref.impulse[impulse_index], 
          exist_state_constraint);
      split_ocps_ref.impulse[impulse_index].computeCondensedPrimalDirection(
          robot_ref, s.impulse[impulse_index], 
          d_ref.impulse[impulse_index]);
      if (split_ocps_ref.impulse[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref.impulse[impulse_index].maxPrimalStepSize();
      if (split_ocps_ref.impulse[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref.impulse[impulse_index].maxDualStepSize();
      // aux 
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.aux[impulse_index], d_ref.aux[impulse_index], 
          exist_state_constraint);
      factorizer.aux[impulse_index].computeControlInputDirection(
          factorization[i+1], d_ref.aux[impulse_index], 
          exist_state_constraint);
      split_ocps_ref.aux[impulse_index].computeCondensedPrimalDirection(
          robot_ref, dtau_aux, s.aux[impulse_index], d_ref.aux[impulse_index]);
      if (split_ocps_ref.aux[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref.aux[impulse_index].maxPrimalStepSize();
      if (split_ocps_ref.aux[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref.aux[impulse_index].maxDualStepSize();
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double lift_time = contact_sequence.liftTime(lift_index);
      const double dtau_lift = lift_time - i * dtau;
      const double dtau_aux = dtau - dtau_lift;
      ASSERT_TRUE(dtau_lift > 0);
      ASSERT_TRUE(dtau_aux > 0);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                 exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization.lift[lift_index], d_ref[i], 
          exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(robot_ref, dtau_lift,
                                                        s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
      // lift
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.lift[lift_index], d_ref.lift[lift_index], 
          exist_state_constraint);
      factorizer.lift[lift_index].computeControlInputDirection(
          factorization[i+1], d_ref.lift[lift_index], 
          exist_state_constraint);
      split_ocps_ref.lift[lift_index].computeCondensedPrimalDirection(
          robot_ref, dtau_aux, s.lift[lift_index], d_ref.lift[lift_index]);
      if (split_ocps_ref.lift[lift_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref.lift[lift_index].maxPrimalStepSize();
      if (split_ocps_ref.lift[lift_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref.lift[lift_index].maxDualStepSize();
    }
    else {
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                 exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization[i+1], d_ref[i], exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(
          robot_ref, dtau, s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
    }
  }
  SplitRiccatiFactorizer::computeCostateDirection(factorization[N], d_ref[N], false);
  if (split_ocps_ref.terminal.maxPrimalStepSize() < primal_step_size_ref) 
    primal_step_size_ref = split_ocps_ref.terminal.maxPrimalStepSize();
  if (split_ocps_ref.terminal.maxDualStepSize() < dual_step_size_ref) 
    dual_step_size_ref = split_ocps_ref.terminal.maxDualStepSize();
  testIsSame(d, d_ref);
  EXPECT_DOUBLE_EQ(primal_step_size, primal_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
}


TEST_F(RiccatiDirectionCalculatorTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  test(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(RiccatiDirectionCalculatorTest, floatingBase) {
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
