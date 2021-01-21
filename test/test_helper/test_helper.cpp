#include "test_helper.hpp"

#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/time_varying_configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
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
  if (robot.max_dimf() > 0) {
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
  else {
    return ContactSequence(robot, N);
  }
}
  

std::shared_ptr<CostFunction> CreateCost(const Robot& robot) {
  auto config_cost = std::make_shared<ConfigurationSpaceCost>(robot);
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
  const Eigen::VectorXd qi_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd vi_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd dvi_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_ref(v_ref);
  config_cost->set_a_weight(a_weight);
  config_cost->set_u_weight(u_weight);
  config_cost->set_u_ref(u_ref);
  config_cost->set_qf_weight(qf_weight);
  config_cost->set_vf_weight(vf_weight);
  config_cost->set_qi_weight(qi_weight);
  config_cost->set_vi_weight(vi_weight);
  config_cost->set_dvi_weight(dvi_weight);
  const int task_frame = 10;
  auto contact_force_cost = std::make_shared<ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, fi_weight;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
    fi_weight.push_back(Eigen::Vector3d::Constant(0.005));
  }
  contact_force_cost->set_f_weight(f_weight);
  contact_force_cost->set_fi_weight(fi_weight);
  auto time_varying_config_cost = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot);
  const double t_begin = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double t_end = t_begin + std::abs(Eigen::VectorXd::Random(1)[0]);
  time_varying_config_cost->set_ref(robot, t_begin, t_end, q_ref, v_ref);
  time_varying_config_cost->set_q_weight(q_weight);
  time_varying_config_cost->set_v_weight(v_weight);
  time_varying_config_cost->set_a_weight(a_weight);
  time_varying_config_cost->set_qf_weight(qf_weight);
  time_varying_config_cost->set_vf_weight(vf_weight);
  time_varying_config_cost->set_qf_weight(qf_weight);
  time_varying_config_cost->set_vf_weight(vf_weight);
  time_varying_config_cost->set_qi_weight(qi_weight);
  time_varying_config_cost->set_vi_weight(vi_weight);
  time_varying_config_cost->set_dvi_weight(dvi_weight);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(config_cost);
  cost->push_back(contact_force_cost);
  cost->push_back(time_varying_config_cost);
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


Solution CreateSolution(const Robot& robot, const int N, const int max_num_impulse) {
  Solution s(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution CreateSolution(const Robot& robot, const ContactSequence& contact_sequence, 
                        const double T, const int N, const int max_num_impulse, 
                        const double t, const bool is_parnmpc) {
  if (robot.maxPointContacts() == 0) {
    return CreateSolution(robot, N, max_num_impulse);
  }
  else if (is_parnmpc) {
    ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_impulse);
    parnmpc_discretizer.discretizeOCP(contact_sequence, t);
    Solution s(robot, N, max_num_impulse);
    for (int i=0; i<N; ++i) {
      s[i].setRandom(robot, contact_sequence.contactStatus(parnmpc_discretizer.contactPhase(i)));
    }
    const int num_impulse = contact_sequence.numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      s.impulse[i].setRandom(robot, contact_sequence.impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      s.aux[i].setRandom(robot, contact_sequence.contactStatus(parnmpc_discretizer.contactPhaseBeforeImpulse(i)));
    }
    const int num_lift = contact_sequence.numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      s.lift[i].setRandom(robot, contact_sequence.contactStatus(parnmpc_discretizer.contactPhaseBeforeLift(i)));
    }
    return s;
  }
  else {
    OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
    ocp_discretizer.discretizeOCP(contact_sequence, t);
    Solution s(robot, N, max_num_impulse);
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


Direction CreateDirection(const Robot& robot, const int N, const int max_num_impulse) {
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    d[i].setRandom();
  }
  return d;
}


Direction CreateDirection(const Robot& robot, const ContactSequence& contact_sequence, 
                          const double T, const int N, const int max_num_impulse, const double t) {
  if (robot.maxPointContacts() == 0) {
    return CreateDirection(robot, N, max_num_impulse);
  }
  else {
    OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
    ocp_discretizer.discretizeOCP(contact_sequence, t);
    Direction d(robot, N, max_num_impulse);
    for (int i=0; i<=N; ++i) {
      d[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)));
    }
    const int num_impulse = contact_sequence.numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      d.impulse[i].setRandom(contact_sequence.impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      d.aux[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterImpulse(i)));
    }
    const int num_lift = contact_sequence.numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      d.lift[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterLift(i)));
    }
    return d;
  }
}


KKTMatrix CreateKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence, 
                          const int N, const int max_num_impulse) {
  KKTMatrix kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  for (int i=0; i<=N; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix[i].Fqq().setIdentity();
      kkt_matrix[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix[i].Fvq().setRandom();
    kkt_matrix[i].Fvv().setRandom();
    kkt_matrix[i].Fvu().setRandom();
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.impulse[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    if (robot.hasFloatingBase()) {
      kkt_matrix.impulse[i].Fqq().setIdentity();
      kkt_matrix.impulse[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.impulse[i].Fvq().setRandom();
    kkt_matrix.impulse[i].Fvv().setRandom();
    kkt_matrix.impulse[i].Pq().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix.aux[i].Fqq().setIdentity();
      kkt_matrix.aux[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.aux[i].Fvq().setRandom();
    kkt_matrix.aux[i].Fvv().setRandom();
    kkt_matrix.aux[i].Fvu().setRandom();
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.hasFloatingBase()) {
      kkt_matrix.lift[i].Fqq().setIdentity();
      kkt_matrix.lift[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.lift[i].Fvq().setRandom();
    kkt_matrix.lift[i].Fvv().setRandom();
    kkt_matrix.lift[i].Fvu().setRandom();
  }
  return kkt_matrix;
}


KKTResidual CreateKKTResidual(const Robot& robot, const ContactSequence& contact_sequence, 
                              const int N, const int max_num_impulse) {
  KKTResidual kkt_residual = KKTResidual(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    kkt_residual[i].lx().setRandom();
    kkt_residual[i].lu().setRandom();
    kkt_residual[i].Fx().setRandom();
  }
  const int num_impulse = contact_sequence.numImpulseEvents();
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    kkt_residual.impulse[i].lx().setRandom();
    kkt_residual.impulse[i].Fx().setRandom();
    kkt_residual.impulse[i].P().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.aux[i].lx().setRandom();
    kkt_residual.aux[i].lu().setRandom();
    kkt_residual.aux[i].Fx().setRandom();
  }
  const int num_lift = contact_sequence.numLiftEvents();
  for (int i=0; i<num_lift; ++i) {
    kkt_residual.lift[i].lx().setRandom();
    kkt_residual.lift[i].lu().setRandom();
    kkt_residual.lift[i].Fx().setRandom();
  }
  return kkt_residual;
}

} // namespace testhelper
} // namespace idocp