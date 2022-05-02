#include "robotoc/mpc/mpc_flying_trotting.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCFlyingTrotting::MPCFlyingTrotting(const Robot& robot, const double T, 
                                     const int N, const int max_steps, 
                                     const int nthreads)
  : foot_step_planner_(),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(robot, max_steps)),
    cost_(std::make_shared<CostFunction>()),
    constraints_(std::make_shared<Constraints>(1.0e-03, 0.995)),
    ocp_solver_(OCP(robot, cost_, constraints_, T, N, max_steps), 
                contact_sequence_, SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(robot.createContactStatus()),
    cs_lfrh_(robot.createContactStatus()),
    cs_rflh_(robot.createContactStatus()),
    cs_flying_(robot.createContactStatus()),
    swing_height_(0),
    flying_time_(0),
    stance_time_(0),
    swing_start_time_(0),
    T_(T),
    dt_(T/N),
    dtm_(T/N),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(N),
    current_step_(0),
    predict_step_(0) {
  try {
    if (robot.maxNumPointContacts() < 4) {
      throw std::out_of_range(
          "invalid argument: input robot is not a quadrupedal robot!\n robot.maxNumPointContacts() must be larger than 4!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  // create costs
  config_cost_ = std::make_shared<ConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.001);
  q_weight.template head<6>().setZero();
  Eigen::VectorXd qi_weight = Eigen::VectorXd::Constant(robot.dimv(), 1);
  qi_weight.template head<6>().setZero();
  config_cost_->set_q_weight(q_weight);
  config_cost_->set_qf_weight(q_weight);
  config_cost_->set_qi_weight(qi_weight);
  config_cost_->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_u_weight(Eigen::VectorXd::Constant(robot.dimu(), 1.0e-02));
  config_cost_->set_vi_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_dvi_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0e-03));
  base_rot_cost_ = std::make_shared<TimeVaryingConfigurationSpaceCost>(robot, base_rot_ref_);
  Eigen::VectorXd base_rot_weight = Eigen::VectorXd::Zero(robot.dimv());
  base_rot_weight.template head<6>() << 0, 0, 0, 1000, 1000, 1000;
  base_rot_cost_->set_q_weight(base_rot_weight);
  base_rot_cost_->set_qf_weight(base_rot_weight);
  base_rot_cost_->set_qi_weight(base_rot_weight);
  LF_foot_cost_ = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, 
                                                               robot.contactFrames()[0],
                                                               LF_foot_ref_);
  LH_foot_cost_ = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, 
                                                               robot.contactFrames()[1],
                                                               LH_foot_ref_);
  RF_foot_cost_ = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, 
                                                               robot.contactFrames()[2],
                                                               RF_foot_ref_);
  RH_foot_cost_ = std::make_shared<TimeVaryingTaskSpace3DCost>(robot, 
                                                               robot.contactFrames()[3],
                                                               RH_foot_ref_);
  LF_foot_cost_->set_x3d_weight(Eigen::Vector3d::Constant(1.0e04));
  LH_foot_cost_->set_x3d_weight(Eigen::Vector3d::Constant(1.0e04));
  RF_foot_cost_->set_x3d_weight(Eigen::Vector3d::Constant(1.0e04));
  RH_foot_cost_->set_x3d_weight(Eigen::Vector3d::Constant(1.0e04));
  com_cost_ = std::make_shared<TimeVaryingCoMCost>(robot, com_ref_);
  com_cost_->set_com_weight(Eigen::Vector3d::Constant(1.0e05));
  cost_->push_back(config_cost_);
  cost_->push_back(base_rot_cost_);
  cost_->push_back(LF_foot_cost_);
  cost_->push_back(LH_foot_cost_);
  cost_->push_back(RF_foot_cost_);
  cost_->push_back(RH_foot_cost_);
  cost_->push_back(com_cost_);
  // create constraints 
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  const double mu = 0.5;
  friction_cone_ = std::make_shared<robotoc::FrictionCone>(robot, mu);
  // auto impulse_friction_cone = std::make_shared<robotoc::ImpulseFrictionCone>(robot, mu);
  constraints_->push_back(joint_position_lower);
  constraints_->push_back(joint_position_upper);
  constraints_->push_back(joint_velocity_lower);
  constraints_->push_back(joint_velocity_upper);
  constraints_->push_back(joint_torques_lower);
  constraints_->push_back(joint_torques_upper);
  constraints_->push_back(friction_cone_);
  // constraints_->push_back(impulse_friction_cone);
  // create contact status
  cs_standing_.activateContacts({0, 1, 2, 3});
  cs_lfrh_.activateContacts({0, 3});
  cs_rflh_.activateContacts({1, 2});
}


MPCFlyingTrotting::MPCFlyingTrotting() {
}


MPCFlyingTrotting::~MPCFlyingTrotting() {
}


void MPCFlyingTrotting::setGaitPattern(
    const std::shared_ptr<FootStepPlannerBase>& foot_step_planner, 
    const double swing_height, const double flying_time, 
    const double stance_time, const double swing_start_time) {
  try {
    if (swing_height <= 0) {
      throw std::out_of_range("invalid value: swing_height must be positive!");
    }
    if (flying_time <= 0) {
      throw std::out_of_range("invalid value: flying_time must be positive!");
    }
    if (stance_time < 0) {
      throw std::out_of_range("invalid value: stance_time must be non-negative!");
    }
    if (swing_start_time <= 0) {
      throw std::out_of_range("invalid value: swing_start_time must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  foot_step_planner_ = foot_step_planner;
  flying_time_ = flying_time;
  stance_time_ = stance_time;
  swing_start_time_ = swing_start_time;
  const int num_phases_in_period = 2;
  LF_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(0, swing_height, 
                                                           swing_start_time_+stance_time_, 
                                                           2.*flying_time_+stance_time_, 
                                                           stance_time_, num_phases_in_period);
  LH_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(1, swing_height, 
                                                           swing_start_time_-flying_time_, 
                                                           2.*flying_time_+stance_time_, 
                                                           stance_time_, num_phases_in_period);
  RF_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(2, swing_height, 
                                                           swing_start_time_-flying_time_, 
                                                           2.*flying_time_+stance_time_, 
                                                           stance_time_, num_phases_in_period);
  RH_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(3, swing_height, 
                                                           swing_start_time_+stance_time_, 
                                                           2.*flying_time_+stance_time_, 
                                                           stance_time_, num_phases_in_period);
  LF_foot_cost_->set_x3d_ref(LF_foot_ref_);
  LH_foot_cost_->set_x3d_ref(LH_foot_ref_);
  RF_foot_cost_->set_x3d_ref(RF_foot_ref_);
  RH_foot_cost_->set_x3d_ref(RH_foot_ref_);
  com_ref_ = std::make_shared<MPCPeriodicCoMRef>(swing_start_time_, 
                                                 flying_time_+stance_time_, 0.0,
                                                 num_phases_in_period);
  com_cost_->set_com_ref(com_ref_);
}


void MPCFlyingTrotting::init(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, 
                             const SolverOptions& solver_options) {
  try {
    if (t >= swing_start_time_) {
      throw std::out_of_range(
          "invalid value: t must be less than" + std::to_string(swing_start_time_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  current_step_ = 0;
  predict_step_ = 0;
  contact_sequence_->initContactSequence(cs_standing_);
  bool add_step = addStep(t);
  while (add_step) {
    add_step = addStep(t);
  }
  foot_step_planner_->init(q);
  config_cost_->set_q_ref(q);
  const int num_phases_in_period = 2;
  base_rot_ref_ = std::make_shared<MPCPeriodicConfigurationRef>(q, swing_start_time_, 
                                                                flying_time_+stance_time_, 0.0,
                                                                num_phases_in_period);
  base_rot_cost_->set_q_ref(base_rot_ref_);
  resetContactPlacements(q, v);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  ts_last_ = swing_start_time_;
}


void MPCFlyingTrotting::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCFlyingTrotting::updateSolution(const double t, const double dt,
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  assert(dt > 0);
  const bool add_step = addStep(t);
  const auto ts = contact_sequence_->eventTimes();
  bool remove_step = false;
  if (!ts.empty()) {
    if (ts.front()+eps_ < t+dt) {
      ts_last_ = ts.front();
      ocp_solver_.extrapolateSolutionInitialPhase(t);
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetContactPlacements(q, v);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCFlyingTrotting::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


const Solution& MPCFlyingTrotting::getSolution() const {
  return ocp_solver_.getSolution();
}


const hybrid_container<LQRPolicy>& MPCFlyingTrotting::getLQRPolicy() const {
  return ocp_solver_.getLQRPolicy();
}


double MPCFlyingTrotting::KKTError(const double t, const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCFlyingTrotting::KKTError() const {
  return ocp_solver_.KKTError();
}


std::shared_ptr<CostFunction> MPCFlyingTrotting::getCostHandle() {
  return cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCFlyingTrotting::getConfigCostHandle() {
  return config_cost_;
}


std::shared_ptr<TimeVaryingConfigurationSpaceCost> MPCFlyingTrotting::getBaseRotationCostHandle() {
  return base_rot_cost_;
}


std::vector<std::shared_ptr<TimeVaryingTaskSpace3DCost>> MPCFlyingTrotting::getSwingFootCostHandle() {
  std::vector<std::shared_ptr<TimeVaryingTaskSpace3DCost>> swing_foot_cost;
  swing_foot_cost = {LF_foot_cost_, LH_foot_cost_, RF_foot_cost_, RH_foot_cost_};
  return swing_foot_cost;
}


std::shared_ptr<TimeVaryingCoMCost> MPCFlyingTrotting::getCoMCostHandle() {
  return com_cost_;
}


std::shared_ptr<Constraints> MPCFlyingTrotting::getConstraintsHandle() {
  return constraints_;
}


std::shared_ptr<FrictionCone> MPCFlyingTrotting::getFrictionConeHandle() {
  return friction_cone_;
}


std::shared_ptr<ContactSequence> MPCFlyingTrotting::getContactSequenceHandle() {
  return contact_sequence_;
}


bool MPCFlyingTrotting::addStep(const double t) {
  if (predict_step_ == 0) {
    if (swing_start_time_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_lfrh_, swing_start_time_);
      ++predict_step_;
      return true;
    }
  }
  else {
    double tt = ts_last_;
    if (current_step_%2 == 0) {
      tt += flying_time_;
    }
    else {
      tt += stance_time_;
    }
    const auto ts = contact_sequence_->eventTimes();
    if (!ts.empty()) {
      if (predict_step_%2 == 0) {
        tt = ts.back() + flying_time_;
      }
      else {
        tt = ts.back() + stance_time_;
      }
    }
    if (tt < t+T_-dtm_) {
      if (predict_step_%4 == 0) {
        contact_sequence_->push_back(cs_lfrh_, tt);
      }
      else if (predict_step_%4 == 2) {
        contact_sequence_->push_back(cs_rflh_, tt);
      }
      else {
        contact_sequence_->push_back(cs_flying_, tt);
      }
      ++predict_step_;
      return true;
    }
  }
  return false;
}


void MPCFlyingTrotting::resetContactPlacements(const Eigen::VectorXd& q,
                                               const Eigen::VectorXd& v) {
  const bool success = foot_step_planner_->plan(q, v, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases()+1);
  for (int phase=0; phase<contact_sequence_->numContactPhases(); ++phase) {
    contact_sequence_->setContactPlacements(phase, 
                                            foot_step_planner_->contactPosition(phase+1));
  }
  base_rot_ref_->setConfigurationRef(contact_sequence_, foot_step_planner_);
  LF_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  LH_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  RF_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  RH_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  com_ref_->setCoMRef(contact_sequence_, foot_step_planner_);
}

} // namespace robotoc 