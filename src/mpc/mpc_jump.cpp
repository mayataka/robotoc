#include "robotoc/mpc/mpc_jump.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>


namespace robotoc {

MPCJump::MPCJump(const Robot& robot, const double T, const int N, 
                 const int nthreads)
  : foot_step_planner_(),
    contact_sequence_(std::make_shared<robotoc::ContactSequence>(robot, 1)),
    cost_(std::make_shared<CostFunction>()),
    constraints_(std::make_shared<Constraints>(1.0e-03, 0.995)),
    sto_cost_(std::make_shared<STOCostFunction>()),
    sto_constraints_(std::make_shared<STOConstraints>(2)),
    robot_(robot),
    ocp_solver_(OCP(robot, cost_, constraints_, sto_cost_, sto_constraints_, 
                    contact_sequence_, T, N), 
                SolverOptions::defaultOptions(), nthreads), 
    solver_options_(SolverOptions::defaultOptions()),
    cs_ground_(robot.createContactStatus()),
    cs_flying_(robot.createContactStatus()),
    s_(),
    flying_time_(0),
    min_flying_time_(0),
    ground_time_(0),
    min_ground_time_(0),
    T_(T),
    dt_(T/N),
    dtm_(T/N),
    t_mpc_start_(0),
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    N_(N),
    current_step_(0),
    nthreads_(nthreads) {
  // create costs
  config_cost_ = std::make_shared<ConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_weight = Eigen::VectorXd::Constant(robot.dimv(), 0.01);
  q_weight.template head<6>() << 0, 100, 100, 100, 100, 100;
  Eigen::VectorXd q_weight_impulse = Eigen::VectorXd::Constant(robot.dimv(), 10);
  q_weight_impulse.template head<6>() << 0, 1000, 1000, 1000, 1000, 1000;
  config_cost_->set_q_weight(q_weight);
  config_cost_->set_q_weight_terminal(q_weight);
  config_cost_->set_q_weight_impulse(q_weight_impulse);
  config_cost_->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 1.0));
  config_cost_->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 1.0e-03));
  config_cost_->set_v_weight_impulse(Eigen::VectorXd::Constant(robot.dimv(), 10.0));
  config_cost_->set_dv_weight_impulse(Eigen::VectorXd::Constant(robot.dimv(), 1.0e-03));
  cost_->push_back(config_cost_);
  // create constraints 
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  const double mu = 0.5;
  friction_cone_ = std::make_shared<robotoc::FrictionCone>(robot, mu);
  constraints_->push_back(joint_position_lower);
  constraints_->push_back(joint_position_upper);
  constraints_->push_back(joint_velocity_lower);
  constraints_->push_back(joint_velocity_upper);
  constraints_->push_back(joint_torques_lower);
  constraints_->push_back(joint_torques_upper);
  constraints_->push_back(friction_cone_);
  // init contact status
  for (int i=0; i<cs_ground_.maxNumContacts(); ++i) {
    cs_ground_.activateContact(i);
  }
}


MPCJump::MPCJump() {
}


MPCJump::~MPCJump() {
}


MPCJump MPCJump::clone() const {
  auto mpc = MPCJump(robot_, T_, N_, nthreads_);
  auto contact_planner = foot_step_planner_->clone();
  mpc.setJumpPattern(contact_planner, flying_time_, min_flying_time_, 
                     ground_time_, min_ground_time_);
  mpc.setSolverOptions(solver_options_);
  return mpc;
}


void MPCJump::setJumpPattern(
    const std::shared_ptr<ContactPlannerBase>& foot_step_planner, 
    const double flying_time, const double min_flying_time, 
    const double ground_time, const double min_ground_time) {
  try {
    if (flying_time <= 0) {
      throw std::out_of_range("invalid value: flying_time must be positive!");
    }
    if (min_flying_time <= 0) {
      throw std::out_of_range("invalid value: min_flying_time must be positive!");
    }
    if (ground_time <= 0) {
      throw std::out_of_range("invalid value: ground_time must be positive!");
    }
    if (min_ground_time <= 0) {
      throw std::out_of_range("invalid value: min_ground_time must be positive!");
    }
    if (flying_time+ground_time > T_) {
      throw std::out_of_range("invalid value: flying_time+ground_time must be less than T!");
    }
    if (min_flying_time+min_ground_time > T_) {
      throw std::out_of_range("invalid value: min_flying_time+min_ground_time must be less than T!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  foot_step_planner_ = foot_step_planner;
  flying_time_ = flying_time;
  min_flying_time_ = min_flying_time;
  ground_time_ = ground_time;
  min_ground_time_ = min_ground_time;
}


void MPCJump::init(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, 
                   const SolverOptions& solver_options, const bool sto) {
  current_step_ = 0;
  contact_sequence_->init(cs_ground_);
  const double t_lift_off   = t + T_ - ground_time_ - flying_time_;
  const double t_touch_down = t + T_ - ground_time_;
  contact_sequence_->reserve(1);
  contact_sequence_->push_back(cs_flying_, t_lift_off, sto);
  contact_sequence_->push_back(cs_ground_, t_touch_down, sto);
  resetMinimumDwellTimes(t, dtm_);
  foot_step_planner_->init(q);
  Eigen::VectorXd q_ref = q;
  q_ref.template head<3>().noalias() += (foot_step_planner_->CoM(2) 
                                          - foot_step_planner_->CoM(0));
  q_ref.template segment<4>(3) = Eigen::Quaterniond(foot_step_planner_->R(2)).coeffs();
  config_cost_->set_q_ref(q_ref);
  resetContactPlacements(t, q, v);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.solve(t, q, v, true);
  s_ = ocp_solver_.getSolution();
  const auto ts = contact_sequence_->eventTimes();
  ground_time_ = t + T_ - ts[1];
  flying_time_ = t + T_ - ts[0] - ground_time_;
  t_mpc_start_ = t;
}


void MPCJump::reset(const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, 
                    const SolverOptions& solver_options, const bool sto) {
  current_step_ = 0;
  contact_sequence_->init(cs_ground_);
  const double t_lift_off   = t + T_ - ground_time_ - flying_time_;
  const double t_touch_down = t + T_ - ground_time_;
  contact_sequence_->push_back(cs_flying_, t_lift_off, sto);
  contact_sequence_->push_back(cs_ground_, t_touch_down, sto);
  resetMinimumDwellTimes(t, dtm_);
  foot_step_planner_->init(q);
  Eigen::VectorXd q_ref = q;
  q_ref.template head<3>().noalias() += (foot_step_planner_->CoM(2) 
                                          - foot_step_planner_->CoM(0));
  q_ref.template segment<4>(3) = Eigen::Quaterniond(foot_step_planner_->R(2)).coeffs();
  config_cost_->set_q_ref(q_ref);
  resetContactPlacements(t, q, v);
  ocp_solver_.setSolution(s_);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
  ocp_solver_.setSolverOptions(solver_options);
  ocp_solver_.meshRefinement(t);
  ocp_solver_.solve(t, q, v, true);
  s_ = ocp_solver_.getSolution();
  const auto ts = contact_sequence_->eventTimes();
  ground_time_ = t + T_ - ts[1];
  flying_time_ = t + T_ - ts[0] - ground_time_;
  t_mpc_start_ = t;
}


void MPCJump::reset() {
  ocp_solver_.setSolution(s_);
}


void MPCJump::reset(const Eigen::VectorXd& q, const Eigen::VectorXd& v) {
  ocp_solver_.setSolution(s_);
  ocp_solver_.setSolution("q", q);
  ocp_solver_.setSolution("v", v);
}


void MPCJump::setSolverOptions(const SolverOptions& solver_options) {
  ocp_solver_.setSolverOptions(solver_options);
}


void MPCJump::updateSolution(const double t, const double dt,
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(dt > 0);
  const auto ts = contact_sequence_->eventTimes();
  bool remove_step = false;
  if (!ts.empty()) {
    if (ts.front()+eps_ < t+dt) {
      ocp_solver_.extrapolateSolutionInitialPhase(t);
      contact_sequence_->pop_front();
      remove_step = true;
      ++current_step_;
    }
  }
  resetMinimumDwellTimes(t, dt);
  resetContactPlacements(t, q, v);
  ocp_solver_.solve(t, q, v, true);
}


const Eigen::VectorXd& MPCJump::getInitialControlInput() const {
  return ocp_solver_.getSolution(0).u;
}


const Solution& MPCJump::getSolution() const {
  return ocp_solver_.getSolution();
}


const hybrid_container<LQRPolicy>& MPCJump::getLQRPolicy() const {
  return ocp_solver_.getLQRPolicy();
}


double MPCJump::KKTError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v) {
  return ocp_solver_.KKTError(t, q, v);
}


double MPCJump::KKTError() const {
  return ocp_solver_.KKTError();
}


std::shared_ptr<CostFunction> MPCJump::getCostHandle() {
  return cost_;
}


std::shared_ptr<ConfigurationSpaceCost> MPCJump::getConfigCostHandle() {
  return config_cost_;
}


std::shared_ptr<Constraints> MPCJump::getConstraintsHandle() {
  return constraints_;
}


std::shared_ptr<FrictionCone> MPCJump::getFrictionConeHandle() {
  return friction_cone_;
}


void MPCJump::setRobotProperties(const RobotProperties& properties) {
  ocp_solver_.setRobotProperties(properties);
}


void MPCJump::resetMinimumDwellTimes(const double t, const double min_dt) {
  const int num_switches = contact_sequence_->numDiscreteEvents();
  if (num_switches > 0) {
    std::vector<double> minimum_dwell_times;
    switch (current_step_) {
      case 0:
        minimum_dwell_times.push_back(min_dt);
        minimum_dwell_times.push_back(min_flying_time_);
        break;
      case 1:
        minimum_dwell_times.push_back(min_dt);
        break;
      default:
        // if current_step_ >= 2, num_switches == 0.
        break;
    }
    minimum_dwell_times.push_back(min_ground_time_+(t-t_mpc_start_));
    sto_constraints_->setMinimumDwellTimes(minimum_dwell_times);
  }
}


void MPCJump::resetContactPlacements(const double t, const Eigen::VectorXd& q,
                                     const Eigen::VectorXd& v) {
  const bool success = foot_step_planner_->plan(t, q, v, contact_sequence_->contactStatus(0),
                                                contact_sequence_->numContactPhases());
  if (current_step_ == 0) {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacements(1));
    contact_sequence_->setContactPlacements(1, foot_step_planner_->contactPlacements(2));
    contact_sequence_->setContactPlacements(2, foot_step_planner_->contactPlacements(2));
  }
  else if (current_step_ == 1) {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacements(1));
    contact_sequence_->setContactPlacements(1, foot_step_planner_->contactPlacements(1));
  }
  else {
    contact_sequence_->setContactPlacements(0, foot_step_planner_->contactPlacements(1));
  }
}

} // namespace robotoc 