#ifndef IDOCP_TEST_HELPER_HPP_
#define IDOCP_TEST_HELPER_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/discrete_event.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/time_varying_configuration_cost.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

namespace idocp {
namespace testhelper {

ContactSequence CreateContactSequence(const Robot& robot, const int N, 
                                      const int max_num_impulse,
                                      const double t0,
                                      const double event_period) {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, N);
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
    tmp.eventTime = t0 + i * event_period + 0.1 * event_period * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.pushBackDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}
  

std::shared_ptr<CostFunction> CreateCost(const Robot& robot) {
  auto joint_cost = std::make_shared<JointSpaceCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd q_ref = robot.generateFeasibleConfiguration();
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
  auto contact_force_cost = std::make_shared<ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
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
  auto impulse_force_cost = std::make_shared<ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  auto time_varying_cost = std::make_shared<TimeVaryingConfigurationCost>(robot);
  const double t0 = Eigen::VectorXd::Random(1)[0];
  time_varying_cost->set_ref(t0, q_ref, v_ref);
  time_varying_cost->set_q_weight(q_weight);
  time_varying_cost->set_v_weight(v_weight);
  time_varying_cost->set_a_weight(a_weight);
  time_varying_cost->set_qf_weight(qf_weight);
  time_varying_cost->set_vf_weight(vf_weight);
  auto impulse_time_varying_cost = std::make_shared<ImpulseTimeVaryingConfigurationCost>(robot);
  impulse_time_varying_cost->set_ref(t0, q_ref, v_ref);
  impulse_time_varying_cost->set_q_weight(q_weight);
  impulse_time_varying_cost->set_v_weight(v_weight);
  impulse_time_varying_cost->set_dv_weight(a_weight);
  cost->push_back(joint_cost);
  cost->push_back(contact_force_cost);
  cost->push_back(impulse_joint_cost);
  cost->push_back(impulse_force_cost);
  cost->push_back(time_varying_cost);
  cost->push_back(impulse_time_varying_cost);
  return cost;
}


std::shared_ptr<Constraints> CreateConstraints(const Robot& robot) {
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


Solution CreateSolution(const Robot& robot, const int N, const int max_num_impulse=0) {
  Solution s(N, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution CreateSolution(const Robot& robot, const ContactSequence& contact_sequence, 
                        const double T, const int N, const int max_num_impulse, const double t) {
  if (robot.maxPointContacts() == 0) {
    return CreateSolution(robot, N, max_num_impulse);
  }
  else {
    OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
    ocp_discretizer.discretizeOCP(contact_sequence, t);
    Solution s(N, max_num_impulse, robot);
    for (int i=0; i<=N; ++i) {
      s[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)));
    }
    const int num_impulse = contact_sequence.numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      s.impulse[i].setRandom(robot, contact_sequence.impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      s.aux[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterImpulse(i)));
    }
    const int num_lift = contact_sequence.numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      s.lift[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterLift(i)));
    }
    return s;
  }
}


template <typename Type, typename ImpulseType>
bool IsApprox(const hybrid_container<Type, ImpulseType>& rhs, 
              const hybrid_container<Type, ImpulseType>& lhs) {
  assert(rhs.data.size() == lhs.data.size());
  assert(rhs.impulse.size() == lhs.impulse.size());
  assert(rhs.aux.size() == lhs.aux.size());
  assert(rhs.lift.size() == lhs.lift.size());
  const int N = rhs.data.size()-1;
  const int max_num_impulse = rhs.impulse.size();
  for (int i=0; i<=N; ++i) {
    if (!rhs[i].isApprox(lhs[i])) {
      return false;
    } 
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.impulse[i].isApprox(lhs.impulse[i])) {
      return false;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.aux[i].isApprox(lhs.aux[i])) {
      return false;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.lift[i].isApprox(lhs.lift[i])) {
      return false;
    }
  }
  return true;
}


template <typename Type, typename ImpulseType>
bool HasNaN(const hybrid_container<Type, ImpulseType>& obj) {
  const int N = obj.data.size()-1;
  const int max_num_impulse = obj.impulse.size();
  for (int i=0; i<=N; ++i) {
    if (obj[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.impulse[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.aux[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.lift[i].hasNaN()) {
      return true;
    }
  }
  return false;
}

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_HPP_