#include "robotoc/mpc/mpc_biped_walk.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCBipedWalk::MPCBipedWalk(const Robot& robot, const double T, const int N, 
                           const int nthreads)
  : foot_step_planner_(),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(robot)),
    cost_(std::make_shared<CostFunction>()),
    constraints_(std::make_shared<Constraints>(1.0e-03, 0.995)),
    ocp_solver_(OCP(robot, cost_, constraints_, contact_sequence_, T, N), 
                SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_standing_(robot.createContactStatus()),
    cs_right_swing_(robot.createContactStatus()),
    cs_left_swing_(robot.createContactStatus()),
    step_height_(0),
    swing_time_(0),
    double_support_time_(0),
    swing_start_time_(0),
    T_(T),
    dt_(T/N),
    dtm_(T/N),
    ts_last_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(N),
    current_step_(0),
    predict_step_(0),
    enable_double_support_phase_(false) {
  if (robot.maxNumSurfaceContacts() < 2) {
    throw std::out_of_range(
        "[MPCBipedWalk] invalid argument: 'robot' is not a bipedal robot!\n robot.maxNumSurfaceContacts() must be larger than 2!");
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
  config_cost_->set_dv_weight_impulse(Eigen::VectorXd::Constant(robot.dimv(), 1.0e-02));
  base_rot_cost_ = std::make_shared<ConfigurationSpaceCost>(robot, base_rot_ref_);
  Eigen::VectorXd base_rot_weight = Eigen::VectorXd::Zero(robot.dimv());
  base_rot_weight.template head<6>() << 0, 0, 0, 1000, 1000, 1000;
  base_rot_cost_->set_q_weight(base_rot_weight);
  base_rot_cost_->set_q_weight_terminal(base_rot_weight);
  base_rot_cost_->set_q_weight_impulse(base_rot_weight);
  L_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[0],
                                                   L_foot_ref_);
  R_foot_cost_ = std::make_shared<TaskSpace3DCost>(robot, robot.contactFrames()[1],
                                                   R_foot_ref_);
  L_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  R_foot_cost_->set_weight(Eigen::Vector3d::Constant(1.0e04));
  com_cost_ = std::make_shared<CoMCost>(robot, com_ref_);
  com_cost_->set_weight(Eigen::Vector3d::Constant(1.0e03));
  cost_->push_back(config_cost_);
  cost_->push_back(base_rot_cost_);
  cost_->push_back(L_foot_cost_);
  cost_->push_back(R_foot_cost_);
  cost_->push_back(com_cost_);
  // create constraints 
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  const double X = 0.1;
  const double Y = 0.05;
  contact_wrench_cone_ = std::make_shared<robotoc::ContactWrenchCone>(robot, X, Y);
  impulse_wrench_cone_ = std::make_shared<robotoc::ImpulseWrenchCone>(robot, X, Y);
  constraints_->push_back(joint_position_lower);
  constraints_->push_back(joint_position_upper);
  constraints_->push_back(joint_velocity_lower);
  constraints_->push_back(joint_velocity_upper);
  constraints_->push_back(joint_torques_lower);
  constraints_->push_back(joint_torques_upper);
  constraints_->push_back(contact_wrench_cone_);
  constraints_->push_back(impulse_wrench_cone_);
  // create contact status
  cs_standing_.activateContacts(std::vector<int>({0, 1}));
  cs_right_swing_.activateContacts(std::vector<int>({0}));
  cs_left_swing_.activateContacts(std::vector<int>({1}));
  const double friction_coefficient = 0.5;
  cs_standing_.setFrictionCoefficients(std::vector<double>(2, friction_coefficient));
  cs_right_swing_.setFrictionCoefficients(std::vector<double>(2, friction_coefficient));
  cs_left_swing_.setFrictionCoefficients(std::vector<double>(2, friction_coefficient));
}


MPCBipedWalk::MPCBipedWalk() {
}


MPCBipedWalk::~MPCBipedWalk() {
}


void MPCBipedWalk::setGaitPattern(const std::shared_ptr<ContactPlannerBase>& foot_step_planner,
                                  const double swing_height, const double swing_time,
                                  const double double_support_time,
                                  const double swing_start_time) {
  if (swing_height <= 0) {
    throw std::out_of_range("[MPCBipedWalk] invalid argument: 'swing_height' must be positive!");
  }
  if (swing_time <= 0) {
    throw std::out_of_range("[MPCBipedWalk] invalid argument: 'swing_time' must be positive!");
  }
  if (double_support_time < 0) {
    throw std::out_of_range("[MPCBipedWalk] invalid argument: 'double_support_time' must be non-negative!");
  }
  if (swing_start_time <= 0) {
    throw std::out_of_range("[MPCBipedWalk] invalid argument: 'swing_start_time' must be positive!");
  }
  foot_step_planner_ = foot_step_planner;
  swing_time_ = swing_time;
  double_support_time_ = double_support_time;
  swing_start_time_ = swing_start_time;
  L_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(0, swing_height, 
                                                          swing_start_time_+swing_time_+double_support_time_, 
                                                          swing_time_, swing_time_+2*double_support_time_);
  R_foot_ref_ = std::make_shared<MPCPeriodicSwingFootRef>(1, swing_height, 
                                                          swing_start_time_, 
                                                          swing_time_, swing_time_+2*double_support_time_);
  L_foot_cost_->set_ref(L_foot_ref_);
  R_foot_cost_->set_ref(R_foot_ref_);
  com_ref_ = std::make_shared<MPCPeriodicCoMRef>(swing_start_time_, 
                                                 swing_time_, double_support_time_);
  com_cost_->set_ref(com_ref_);
}


void MPCBipedWalk::init(const double t, const Eigen::VectorXd& q, 
                        const Eigen::VectorXd& v, 
                        const SolverOptions& solver_options) {
  if (t >= swing_start_time_) {
    throw std::out_of_range(
        "[MPCBipedWalk] invalid argument: 't' must be less than " + std::to_string(swing_start_time_) + "!");
  }
  current_step_ = 0;
  predict_step_ = 0;
  contact_sequence_->reserve(std::floor(T_/(swing_time_+double_support_time_)));
  contact_sequence_->init(cs_standing_);
  bool add_step = addStep(t);
  while (add_step) {
    add_step = addStep(t);
  }
  foot_step_planner_->init(q);
  config_cost_->set_q_ref(q);
  base_rot_ref_ = std::make_shared<MPCPeriodicConfigurationRef>(q, swing_start_time_, 
                                                                swing_time_, double_support_time_);
  base_rot_cost_->set_ref(base_rot_ref_);
  resetContactPlacements(t, q, v);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  ts_last_ = swing_start_time_;
}


void MPCBipedWalk::reset() {
  ocp_solver_.setSolution(s_);
}


void MPCBipedWalk::reset(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  ocp_solver_.setSolution(s_);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
}


void MPCBipedWalk::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCBipedWalk::updateSolution(const double t, const double dt,
                                 const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v) {
  assert(dt > 0);
  const bool add_step = addStep(t);
  const auto ts = contact_sequence_->eventTimes();
  bool remove_step = false;
  if (!ts.empty()) {
    if (ts.front()+eps_ < t+dt) {
      ts_last_ = ts.front();
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetContactPlacements(t, q, v);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCBipedWalk::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


const Solution& MPCBipedWalk::getSolution() const {
  return ocp_solver_.getSolution();
}


const hybrid_container<LQRPolicy>& MPCBipedWalk::getLQRPolicy() const {
  return ocp_solver_.getLQRPolicy();
}


double MPCBipedWalk::KKTError(const double t, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCBipedWalk::KKTError() const {
  return ocp_solver_.KKTError();
}


std::shared_ptr<CostFunction> MPCBipedWalk::getCostHandle() {
  return cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCBipedWalk::getConfigCostHandle() {
  return config_cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCBipedWalk::getBaseRotationCostHandle() {
  return base_rot_cost_;
}


std::vector<std::shared_ptr<TaskSpace3DCost>> MPCBipedWalk::getSwingFootCostHandle() {
  std::vector<std::shared_ptr<TaskSpace3DCost>> swing_foot_cost;
  swing_foot_cost = {L_foot_cost_, R_foot_cost_};
  return swing_foot_cost;
}


std::shared_ptr<CoMCost> MPCBipedWalk::getCoMCostHandle() {
  return com_cost_;
}


std::shared_ptr<Constraints> MPCBipedWalk::getConstraintsHandle() {
  return constraints_;
}


std::shared_ptr<ContactWrenchCone> MPCBipedWalk::getContactWrenchConeHandle() {
  return contact_wrench_cone_;
}


std::shared_ptr<ImpulseWrenchCone> MPCBipedWalk::getImpulseWrenchConeHandle() {
  return impulse_wrench_cone_;
}


void MPCBipedWalk::setRobotProperties(const RobotProperties& properties) {
  ocp_solver_.setRobotProperties(properties);
}


bool MPCBipedWalk::addStep(const double t) {
  if (predict_step_ == 0) {
    if (swing_start_time_ < t+T_-dtm_) {
      contact_sequence_->push_back(cs_right_swing_, swing_start_time_);
      ++predict_step_;
      return true;
    }
  }
  else {
    if (enable_double_support_phase_) {
      double tt = ts_last_;
      if (current_step_%2 == 0) {
        tt += double_support_time_;
      }
      else {
        tt += swing_time_;
      }
      const auto ts = contact_sequence_->eventTimes();
      if (!ts.empty()) {
        if (predict_step_%2 == 0) {
          tt = ts.back() + double_support_time_;
        }
        else {
          tt = ts.back() + swing_time_;
        }
      }
      if (tt < t+T_-dtm_) {
        if (predict_step_%4 == 0) {
          contact_sequence_->push_back(cs_right_swing_, tt);
        }
        else if (predict_step_%4 == 2) {
          contact_sequence_->push_back(cs_left_swing_, tt);
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
          contact_sequence_->push_back(cs_right_swing_, tt);
        }
        else {
          contact_sequence_->push_back(cs_left_swing_, tt);
        }
        ++predict_step_;
        return true;
      }
    }
  }
  return false;
}


void MPCBipedWalk::resetContactPlacements(const double t, const Eigen::VectorXd& q,
                                          const Eigen::VectorXd& v) {
  const bool success = foot_step_planner_->plan(t, q, v, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases());
  for (int phase=0; phase<contact_sequence_->numContactPhases(); ++phase) {
    contact_sequence_->setContactPlacements(phase, 
                                            foot_step_planner_->contactPlacements(phase+1));
  }
  base_rot_ref_->setConfigurationRef(contact_sequence_, foot_step_planner_);
  L_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  R_foot_ref_->setSwingFootRef(contact_sequence_, foot_step_planner_);
  com_ref_->setCoMRef(contact_sequence_, foot_step_planner_);
}

} // namespace robotoc 