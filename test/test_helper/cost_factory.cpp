#include "cost_factory.hpp"

#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/local_contact_force_cost.hpp"


namespace robotoc {
namespace testhelper {

ConfigurationSpaceRef::ConfigurationSpaceRef(
    const Eigen::VectorXd& q0_ref, const Eigen::VectorXd& v_ref)
  : q0_ref_(q0_ref),
    v_ref_(v_ref) {
}


void ConfigurationSpaceRef::updateRef(const Robot& robot, 
                                      const GridInfo& grid_info,
                                      Eigen::VectorXd& q_ref) const {
  robot.integrateConfiguration(q0_ref_, v_ref_, grid_info.t, q_ref);
}


bool ConfigurationSpaceRef::isActive(const GridInfo& grid_info) const {
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
  const Eigen::VectorXd q_weight_terminal = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_weight_terminal = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd q_weight_impact = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_weight_impact = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd dv_weight_impact = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_ref(v_ref);
  config_cost->set_a_weight(a_weight);
  config_cost->set_u_weight(u_weight);
  config_cost->set_u_ref(u_ref);
  config_cost->set_q_weight_terminal(q_weight_terminal);
  config_cost->set_v_weight_terminal(v_weight_terminal);
  config_cost->set_q_weight_impact(q_weight_impact);
  config_cost->set_v_weight_impact(v_weight_impact);
  config_cost->set_dv_weight_impact(dv_weight_impact);
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
  auto time_varying_config_ref = std::make_shared<ConfigurationSpaceRef>(q0_ref, v0_ref);
  auto time_varying_config_cost 
      = std::make_shared<ConfigurationSpaceCost>(robot, time_varying_config_ref);
  time_varying_config_cost->set_q_weight(q_weight);
  time_varying_config_cost->set_q_weight_terminal(q_weight_terminal);
  time_varying_config_cost->set_q_weight_impact(q_weight_impact);
  cost->push_back(time_varying_config_cost);

  return cost;
}

} // namespace testhelper
} // namespace robotoc