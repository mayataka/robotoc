#include "robotoc/mpc/mpc_trot.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCTrot::MPCTrot(const Robot& robot, const double T, const int N, 
                 const int nthreads)
  : foot_step_planner_(),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(robot)),
    cost_(std::make_shared<CostFunction>()),
    constraints_(std::make_shared<Constraints>(1.0e-03, 0.995)),
    ocp_solver_(OCP(robot, cost_, constraints_, contact_sequence_, T, N), 
                SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(robot.createContactStatus()),
    cs_lfrh_(robot.createContactStatus()),
    cs_rflh_(robot.createContactStatus()),
    swing_height_(0),
    swing_time_(0),
    stance_time_(0),
    swing_start_time_(0),
    T_(T),
    dt_(T/N),
    dtm_(T/N),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(N),
    current_step_(0),
    predict_step_(0),
    enable_stance_phase_(false) {
  if (robot.maxNumPointContacts() < 4) {
    throw std::out_of_range(
        "[MPCTrot] invalid argument: 'robot' is not a quadrupedal robot!\n robot.maxNumPointContacts() must be larger than 4!");
  }
  // create costs
  config_cost_ = std::make_shared<ConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.001);
  q_weight.template head<6>().setZero();
  Eigen::VectorXd q_weight_impulse = Eigen::VectorXd::Constant(robot.dimv(), 1);
  q_weight_impulse.template head<6>().setZero();
  config_cost_->set_q_weight(q_weight);
  config_cost_->set_q_weight_terminal(q_weight);
  config_cost_->set_q_weight_impulse(q_weight_impulse);
  config_cost_->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_u_weight(Eigen::VectorXd::Constant(robot.dimu(), 1.0e-02));
  config_cost_->set_v_weight_impulse(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_dv_weight_impulse(Eigen::VectorXd::Constant(robot.dimv(), 1.0e-03));
  base_rot_cost_ = std::make_shared<ConfigurationSpaceCost>(robot, base_rot_ref_);
  Eigen::VectorXd base_rot_weight = Eigen::VectorXd::Zero(robot.dimv());
  base_rot_weight.template head<6>() << 0, 0, 0, 1000, 1000, 1000;
  base_rot_cost_->set_q_weight(base_rot_weight);
  base_rot_cost_->set_q_weight_terminal(base_rot_weight);
  base_rot_cost_->set_q_weight_impulse(base_rot_weight);
  LF_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[0],
                                                    LF_foot_ref_);
  LH_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[1],
                                                    LH_foot_ref_);
  RF_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[2],
                                                    RF_foot_ref_);
  RH_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[3],
                                                    RH_foot_ref_);
  LF_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  LH_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  RF_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  RH_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  com_cost_ = std::make_shared<CoMCost>(robot, com_ref_);
  com_cost_->set_weight(Eigen::Vector3d::Constant(1.0e03));
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
  cs_standing_.activateContacts(std::vector<int>({0, 1, 2, 3}));
  cs_lfrh_.activateContacts(std::vector<int>({0, 3}));
  cs_rflh_.activateContacts(std::vector<int>({1, 2}));
}


MPCTrot::MPCTrot() {
}


MPCTrot::~MPCTrot() {
}


void MPCTrot::setGaitPattern(const std::shared_ptr<ContactPlannerBase>& foot_step_planner,
                             const double swing_height, const double swing_time,
                             const double stance_time, const double swing_start_time) {
  if (swing_height <= 0) {
    throw std::out_of_range("[MPCTrot] invalid argument: 'swing_height' must be positive!");
  }
  if (swing_time <= 0) {
    throw std::out_of_range("[MPCTrot] invalid argument: 'swing_time' must be positive!");
  }
  if (stance_time < 0) {
    throw std::out_of_range("[MPCTrot] invalid argument: 'stance_time' must be non-negative!");
  }
  if (swing_start_time <= 0) {
    throw std::out_of_range("[MPCTrot] invalid argument: 'swing_start_time' must be positive!");
  }
  foot_step_planner_ = foot_step_planner;
  swing_time_ = swing_time;
  stance_time_ = stance_time;
  swing_start_time_ = swing_start_time;
  enable_stance_phase_ = (stance_time_ > 0.);
  LF_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(0, swing_height, 
                                                           swing_start_time_+swing_time_+stance_time_, 
                                                           swing_time_, swing_time_+2*stance_time_);
  LH_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(1, swing_height, 
                                                           swing_start_time_, 
                                                           swing_time_, swing_time_+2*stance_time_);
  RF_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(2, swing_height, 
                                                           swing_start_time_, 
                                                           swing_time_, swing_time_+2*stance_time_);
  RH_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(3, swing_height, 
                                                           swing_start_time_+swing_time_+stance_time_, 
                                                           swing_time_, swing_time_+2*stance_time_);
  LF_foot_cost_->set_ref(LF_foot_ref_);
  LH_foot_cost_->set_ref(LH_foot_ref_);
  RF_foot_cost_->set_ref(RF_foot_ref_);
  RH_foot_cost_->set_ref(RH_foot_ref_);
  com_ref_ = std::make_shared<MPCPeriodicCoMRef>(swing_start_time_, 
                                                 swing_time_, stance_time_);
  com_cost_->set_ref(com_ref_);
}


void MPCTrot::init(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, 
                   const SolverOptions& solver_options) {
  if (t >= swing_start_time_) {
    throw std::out_of_range(
        "invalid argument: 't' must be less than " + std::to_string(swing_start_time_) + "!");
  }
  current_step_ = 0;
  predict_step_ = 0;
  contact_sequence_->reserve(std::floor(T_/(swing_time_+stance_time_)));
  contact_sequence_->init(cs_standing_);
  bool add_step = addStep(t);
  while (add_step) {
    add_step = addStep(t);
  }
  foot_step_planner_->init(q);
  config_cost_->set_q_ref(q);
  base_rot_ref_ = std::make_shared<MPCPeriodicConfigurationRef>(q, swing_start_time_, 
                                                                swing_time_, stance_time_);
  base_rot_cost_->set_ref(base_rot_ref_);
  resetContactPlacements(t, q, v);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  s_ = ocp_solver_.getSolution();
  ts_last_ = swing_start_time_;
}


void MPCTrot::reset() {
  ocp_solver_.setSolution(s_);
}


void MPCTrot::reset(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  ocp_solver_.setSolution(s_);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
}


void MPCTrot::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCTrot::updateSolution(const double t, const double dt,
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
  resetContactPlacements(t, q, v);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCTrot::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


const Solution& MPCTrot::getSolution() const {
  return ocp_solver_.getSolution();
}


const hybrid_container<LQRPolicy>& MPCTrot::getLQRPolicy() const {
  return ocp_solver_.getLQRPolicy();
}


double MPCTrot::KKTError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCTrot::KKTError() const {
  return ocp_solver_.KKTError();
}


std::shared_ptr<CostFunction> MPCTrot::getCostHandle() {
  return cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCTrot::getConfigCostHandle() {
  return config_cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCTrot::getBaseRotationCostHandle() {
  return base_rot_cost_;
}


std::vector<std::shared_ptr<TaskSpace3DCost>> MPCTrot::getSwingFootCostHandle() {
  std::vector<std::shared_ptr<TaskSpace3DCost>> swing_foot_cost;
  swing_foot_cost = {LF_foot_cost_, LH_foot_cost_, RF_foot_cost_, RH_foot_cost_};
  return swing_foot_cost;
}


std::shared_ptr<CoMCost> MPCTrot::getCoMCostHandle() {
  return com_cost_;
}


std::shared_ptr<Constraints> MPCTrot::getConstraintsHandle() {
  return constraints_;
}


std::shared_ptr<FrictionCone> MPCTrot::getFrictionConeHandle() {
  return friction_cone_;
}


void MPCTrot::setRobotProperties(const RobotProperties& properties) {
  ocp_solver_.setRobotProperties(properties);
}


bool MPCTrot::addStep(const double t) {
  if (predict_step_ == 0) {
    if (swing_start_time_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_lfrh_, swing_start_time_);
      ++predict_step_;
      return true;
    }
  }
  else {
    if (enable_stance_phase_) {
      double tt = ts_last_;
      if (current_step_%2 == 0) {
        tt += stance_time_;
      }
      else {
        tt += swing_time_;
      }
      const auto ts = contact_sequence_->eventTimes();
      if (!ts.empty()) {
        if (predict_step_%2 == 0) {
          tt = ts.back() + stance_time_;
        }
        else {
          tt = ts.back() + swing_time_;
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
          contact_sequence_->push_back(cs_standing_, tt);
        }
        ++predict_step_;
        return true;
      }
    }
    else {
      double tt = ts_last_ + swing_time_;
      const auto ts = contact_sequence_->eventTimes();
      if (!ts.empty()) {
        tt = ts.back() + swing_time_;
      }
      if (tt < t+T_-dtm_) {
        if (predict_step_%2 == 0) {
          contact_sequence_->push_back(cs_lfrh_, tt);
        }
        else {
          contact_sequence_->push_back(cs_rflh_, tt);
        }
        ++predict_step_;
        return true;
      }
    }
  }
  return false;
}


void MPCTrot::resetContactPlacements(const double t, const Eigen::VectorXd& q,
                                     const Eigen::VectorXd& v) {
  const bool success = foot_step_planner_->plan(t, q, v, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases());
  for (int phase=0; phase<contact_sequence_->numContactPhases(); ++phase) {
    contact_sequence_->setContactPlacements(phase, 
                                            foot_step_planner_->contactPositions(phase+1),
                                            foot_step_planner_->contactSurfaces(phase+1));
  }
  base_rot_ref_->setConfigurationRef(contact_sequence_, foot_step_planner_);
  LF_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  LH_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  RF_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  RH_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  com_ref_->setCoMRef(contact_sequence_, foot_step_planner_);
}

} // namespace robotoc 