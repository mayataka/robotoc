#include "cost_factory.hpp"

#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/time_varying_configuration_space_cost.hpp"
#include "robotoc/cost/local_contact_force_cost.hpp"


namespace robotoc {
namespace testhelper {

TimeVaryingConfigurationRef::TimeVaryingConfigurationRef(
    const Eigen::VectorXd& q0_ref, const Eigen::VectorXd& v_ref)
  : q0_ref_(q0_ref),
    v_ref_(v_ref) {
}


void TimeVaryingConfigurationRef::update_q_ref(const Robot& robot, 
                                               const double t, 
                                               Eigen::VectorXd& q_ref) const {
  robot.integrateConfiguration(q0_ref_, v_ref_, t, q_ref);
}


bool TimeVaryingConfigurationRef::isActive(const double t) const {
  return true;
}


std::shared_ptr<CostFunction> CreateCost(const Robot& robot) {
  auto cost = std::make_shared<CostFunction>();

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
  cost->push_back(config_cost);

  if (robot.maxNumContacts() > 0) {
    auto local_contact_force_cost = std::make_shared<LocalContactForceCost>(robot);
    std::vector<Eigen::Vector3d> f_weight, fi_weight;
    for (int i=0; i<robot.maxNumContacts(); ++i) {
      f_weight.push_back(Eigen::Vector3d::Constant(0.001));
      fi_weight.push_back(Eigen::Vector3d::Constant(0.005));
    }
    local_contact_force_cost->set_f_weight(f_weight);
    local_contact_force_cost->set_fi_weight(fi_weight);
    cost->push_back(local_contact_force_cost);
  }

  const Eigen::VectorXd q0_ref = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v0_ref = Eigen::VectorXd::Random(robot.dimv());
  auto time_varying_config_ref = std::make_shared<TimeVaryingConfigurationRef>(q0_ref, v0_ref);
  auto time_varying_config_cost 
      = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot, time_varying_config_ref);
  time_varying_config_cost->set_q_weight(q_weight);
  time_varying_config_cost->set_qf_weight(qf_weight);
  time_varying_config_cost->set_qi_weight(qi_weight);
  cost->push_back(time_varying_config_cost);

  return cost;
}

} // namespace testhelper
} // namespace robotoc