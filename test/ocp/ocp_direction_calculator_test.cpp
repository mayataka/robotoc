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
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/ocp_direction_calculator.hpp"


namespace idocp {

class OCPDirectionCalculatorTest : public ::testing::Test {
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

  void testComputeDirection(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, t, dtau;
};


std::shared_ptr<CostFunction> OCPDirectionCalculatorTest::createCost(const Robot& robot) {
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


std::shared_ptr<Constraints> OCPDirectionCalculatorTest::createConstraints(const Robot& robot) {
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


HybridSolution OCPDirectionCalculatorTest::createSolution(const Robot& robot) const {
  HybridSolution s(N+1, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


HybridSolution OCPDirectionCalculatorTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
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


ContactSequence OCPDirectionCalculatorTest::createContactSequence(const Robot& robot) const {
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


void OCPDirectionCalculatorTest::testComputeDirection(const Robot& robot) const {
  auto cost = createCost(robot);
  auto constraints = createConstraints(robot);
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
  OCPLinearizer linearizer(T, N, max_num_impulse, nproc);
  std::vector<Robot> robots(nproc, robot);
  linearizer.linearizeOCP(split_ocps, robots, contact_sequence, t, q, v, s, kkt_matrix, kkt_residual);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  HybridRiccatiFactorization riccati_factorization(N, max_num_impulse, robot);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, riccati_factorization);
  riccati_recursion.backwardRiccatiRecursion(contact_sequence, kkt_matrix, kkt_residual, riccati_factorization);
  riccati_recursion.forwardRiccatiRecursionParallel(contact_sequence, kkt_matrix, kkt_residual);
  riccati_recursion.forwardRiccatiRecursionSerial(contact_sequence, kkt_matrix, kkt_residual, riccati_factorization);
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization(max_num_impulse, 
                                                                            StateConstraintRiccatiFactorization(robot, N, max_num_impulse));
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_impulse; ++i) {
    constraint_factorization[i].setImpulseStatus(contact_sequence.impulseStatus(i));
  }
  riccati_recursion.backwardStateConstraintFactorization(contact_sequence, kkt_matrix, constraint_factorization);
  HybridDirection d = HybridDirection(N, max_num_impulse, robot);
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
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  auto split_ocps_ref = split_ocps;
  auto riccati_factorization_ref = riccati_factorization;
  const auto factorizer = riccati_recursion.getFactorizersHandle();
  OCPDirectionCalculator direction_calculator(T, N, max_num_impulse, nproc);
  OCPDirectionCalculator::computeInitialStateDirection(robots, q, v, s, d);
  robot.subtractConfiguration(q, s[0].q, d_ref[0].dq());
  d_ref[0].dv() = v - s[0].v;
  EXPECT_TRUE(d[0].isApprox(d_ref[0]));
  direction_calculator.computeDirection(split_ocps, robots, contact_sequence,  
                                        riccati_recursion.getFactorizersHandle(), 
                                        riccati_factorization, 
                                        constraint_factorization, s, d);
  const double primal_step_size = direction_calculator.maxPrimalStepSize(contact_sequence);
  const double dual_step_size = direction_calculator.maxDualStepSize(contact_sequence);
  auto robot_ref = robot;
  const Eigen::VectorXd dx0 = d_ref[0].dx();
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  double primal_step_size_ref = 1;
  double dual_step_size_ref = 1;
  if (contact_sequence.existImpulseStage(0)) {
    const double dtau_impulse = contact_sequence.impulseTime(0);
    const double dtau_aux = dtau - dtau_impulse;
    ASSERT_TRUE(dtau_impulse > 0);
    ASSERT_TRUE(dtau_impulse < dtau);
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref[0]);
    OCPDirectionCalculator::computePrimalDirectionInitial(
        factorizer[0], riccati_factorization_ref[0], d_ref[0], exist_state_constraint);
    split_ocps_ref[0].computeCondensedPrimalDirection(robot_ref, dtau_impulse, 
                                                      s[0], d_ref[0]);
    if (split_ocps_ref[0].maxPrimalStepSize() < primal_step_size_ref) 
      primal_step_size_ref = split_ocps_ref[0].maxPrimalStepSize();
    if (split_ocps_ref[0].maxDualStepSize() < dual_step_size_ref) 
      dual_step_size_ref = split_ocps_ref[0].maxDualStepSize();
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionImpulse(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref.impulse[0]);
    OCPDirectionCalculator::computePrimalDirectionImpulse(
        riccati_factorization_ref.impulse[0], dx0, d_ref.impulse[0]);
    split_ocps.impulse[0].computeCondensedPrimalDirection(robot_ref, 
                                                          s.impulse[0], 
                                                          d_ref.impulse[0]);
    if (split_ocps_ref.impulse[0].maxPrimalStepSize() < primal_step_size_ref) 
      primal_step_size_ref = split_ocps_ref.impulse[0].maxPrimalStepSize();
    if (split_ocps_ref.impulse[0].maxDualStepSize() < dual_step_size_ref) 
      dual_step_size_ref = split_ocps_ref.impulse[0].maxDualStepSize();
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionAux(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref.aux[0]);
    OCPDirectionCalculator::computePrimalDirection(
        factorizer.aux[0], riccati_factorization_ref.aux[0], dx0, d_ref.aux[0],
        exist_state_constraint);
    if (split_ocps_ref.aux[0].maxPrimalStepSize() < primal_step_size_ref)
      primal_step_size_ref = split_ocps_ref.aux[0].maxPrimalStepSize();
    if (split_ocps_ref.aux[0].maxDualStepSize() < dual_step_size_ref)
      dual_step_size_ref = split_ocps_ref.aux[0].maxDualStepSize();
  }
  else if (contact_sequence.existLiftStage(0)) {
    const double dtau_lift = contact_sequence.liftTime(0);
    const double dtau_aux = dtau - dtau_lift;
    ASSERT_TRUE(dtau_lift > 0);
    ASSERT_TRUE(dtau_lift < dtau);
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref[0]);
    OCPDirectionCalculator::computePrimalDirectionInitial(
        factorizer[0], riccati_factorization_ref[0], d_ref[0], exist_state_constraint);
    split_ocps_ref[0].computeCondensedPrimalDirection(robot_ref, dtau_lift, 
                                                      s[0], d_ref[0]);
    if (split_ocps_ref[0].maxPrimalStepSize() < primal_step_size_ref) 
      primal_step_size_ref = split_ocps_ref[0].maxPrimalStepSize();
    if (split_ocps_ref[0].maxDualStepSize() < dual_step_size_ref) 
      dual_step_size_ref = split_ocps_ref[0].maxDualStepSize();
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionLift(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref.lift[0]);
    OCPDirectionCalculator::computePrimalDirection(
        factorizer.lift[0], riccati_factorization_ref.lift[0], dx0, d_ref.lift[0], 
        exist_state_constraint);
    if (split_ocps_ref.lift[0].maxPrimalStepSize() < primal_step_size_ref) 
      primal_step_size_ref = split_ocps_ref.lift[0].maxPrimalStepSize();
    if (split_ocps_ref.lift[0].maxDualStepSize() < dual_step_size_ref)
      dual_step_size_ref = split_ocps_ref.lift[0].maxDualStepSize();
  }
  else {
    OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
        contact_sequence, constraint_factorization, d_ref.impulse, 0, 
        riccati_factorization_ref[0]);
    OCPDirectionCalculator::computePrimalDirectionInitial(
        factorizer[0], riccati_factorization_ref[0], d_ref[0], exist_state_constraint);
    split_ocps_ref[0].computeCondensedPrimalDirection(robot_ref, dtau, s[0], d_ref[0]);
    if (split_ocps_ref[0].maxPrimalStepSize() < primal_step_size_ref)
      primal_step_size_ref = split_ocps_ref[0].maxPrimalStepSize();
    if (split_ocps_ref[0].maxDualStepSize() < dual_step_size_ref)
      dual_step_size_ref = split_ocps_ref[0].maxDualStepSize();
  }
  for (int i=1; i<N; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double impulse_time = contact_sequence.impulseTime(impulse_index);
      const double dtau_impulse = impulse_time - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
          contact_sequence, constraint_factorization, d_ref.impulse, i, 
          riccati_factorization_ref[i]);
      OCPDirectionCalculator::computePrimalDirection(
          factorizer[i], riccati_factorization_ref[i], dx0, d_ref[i], exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(robot_ref, dtau_impulse, 
                                                        s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionImpulse(
          contact_sequence, constraint_factorization, d_ref.impulse, impulse_index, 
          riccati_factorization_ref.impulse[impulse_index]);
      OCPDirectionCalculator::computePrimalDirectionImpulse(
          riccati_factorization_ref.impulse[impulse_index], dx0, d_ref.impulse[impulse_index]);
      split_ocps.impulse[impulse_index].computeCondensedPrimalDirection(
          robot_ref, s.impulse[impulse_index], d_ref.impulse[impulse_index]);
      if (split_ocps_ref.impulse[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref.impulse[impulse_index].maxPrimalStepSize();
      if (split_ocps_ref.impulse[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref.impulse[impulse_index].maxDualStepSize();
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionAux(
          contact_sequence, constraint_factorization, d_ref.impulse, impulse_index, 
          riccati_factorization_ref.aux[impulse_index]);
      OCPDirectionCalculator::computePrimalDirection(
          factorizer.aux[impulse_index], riccati_factorization_ref.aux[impulse_index], dx0, 
          d_ref.aux[impulse_index], exist_state_constraint);
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
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
          contact_sequence, constraint_factorization, d_ref.impulse, i, 
          riccati_factorization_ref[i]);
      OCPDirectionCalculator::computePrimalDirection(
          factorizer[i], riccati_factorization_ref[i], dx0, d_ref[i], exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(robot_ref, dtau_lift, 
                                                        s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirectionLift(
          contact_sequence, constraint_factorization, d_ref.impulse, lift_index, 
          riccati_factorization_ref.lift[lift_index]);
      OCPDirectionCalculator::computePrimalDirection(
          factorizer.lift[lift_index], riccati_factorization_ref.lift[lift_index], dx0, 
          d_ref.lift[lift_index], exist_state_constraint);
      if (split_ocps_ref.lift[lift_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = split_ocps_ref.lift[lift_index].maxPrimalStepSize();
      if (split_ocps_ref.lift[lift_index].maxDualStepSize() < dual_step_size_ref)
        dual_step_size_ref = split_ocps_ref.lift[lift_index].maxDualStepSize();
    }
    else {
      OCPDirectionCalculator::aggregateLagrangeMultiplierDirection(
          contact_sequence, constraint_factorization, d_ref.impulse, i, 
          riccati_factorization_ref[i]);
      OCPDirectionCalculator::computePrimalDirection(
          factorizer[i], riccati_factorization_ref[i], dx0, d_ref[i], exist_state_constraint);
      split_ocps_ref[i].computeCondensedPrimalDirection(robot_ref, dtau, s[i], d_ref[i]);
      if (split_ocps_ref[i].maxPrimalStepSize() < primal_step_size_ref)
        primal_step_size_ref = split_ocps_ref[i].maxPrimalStepSize();
      if (split_ocps_ref[i].maxDualStepSize() < dual_step_size_ref)
        dual_step_size_ref = split_ocps_ref[i].maxDualStepSize();
    }
  }
  OCPDirectionCalculator::computePrimalDirectionTerminal(riccati_factorization_ref[N], dx0, d_ref[N]);
  if (split_ocps_ref.terminal.maxPrimalStepSize() < primal_step_size_ref)
    primal_step_size_ref = split_ocps_ref.terminal.maxPrimalStepSize();
  if (split_ocps_ref.terminal.maxDualStepSize() < dual_step_size_ref)
    dual_step_size_ref = split_ocps_ref.terminal.maxDualStepSize();
  for (int i=0; i<=N; ++i) {
    d[i].isApprox(d_ref[i]);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.aux[i].isApprox(d_ref.aux[i]);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].isApprox(d_ref.impulse[i]);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.lift[i].isApprox(d_ref.lift[i]);
  }
  EXPECT_DOUBLE_EQ(primal_step_size, primal_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
}


TEST_F(OCPDirectionCalculatorTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testComputeDirection(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  testComputeDirection(robot);
}


TEST_F(OCPDirectionCalculatorTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testComputeDirection(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  testComputeDirection(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
