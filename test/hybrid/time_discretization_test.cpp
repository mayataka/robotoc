#include <random>
#include <vector>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/time_discretization.hpp"
#include "robotoc/robot/robot.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class TimeDiscretizationTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_events = 5;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    T = 1;
    dt = T / N;
    min_dt = std::sqrt(std::numeric_limits<double>::epsilon());
  }

  virtual void TearDown() {
  }

  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;
  std::shared_ptr<ContactSequence> createContactSequenceOnGrid(const Robot& robot) const;

  int N, max_num_events;
  double t, T, dt, min_dt;
};


std::shared_ptr<ContactSequence> TimeDiscretizationTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  auto contact_sequence = std::make_shared<ContactSequence>(robot, max_num_events);
  contact_sequence->initContactSequence(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot.contactTypes());
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_time = t + i * event_period + dt * std::abs(Eigen::VectorXd::Random(1)[0]);
    contact_sequence->push_back(tmp, event_time, false);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


std::shared_ptr<ContactSequence> TimeDiscretizationTest::createContactSequenceOnGrid(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  auto contact_sequence = std::make_shared<ContactSequence>(robot, max_num_events);
  contact_sequence->initContactSequence(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot.contactTypes());
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_time = t + (i+1) * event_period + min_dt * Eigen::VectorXd::Random(1)[0];
    contact_sequence->push_back(tmp, event_time, false);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


TEST_P(TimeDiscretizationTest, constructor) {
  TimeDiscretization discretization(T, N, max_num_events);
  EXPECT_EQ(discretization.N(), N);
  EXPECT_EQ(discretization.N_impulse(), 0);
  EXPECT_EQ(discretization.N_lift(), 0);
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(discretization.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(discretization.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(discretization.isTimeStageBeforeLift(i));
  }
  for (int i=1; i<=N; ++i) {
    EXPECT_FALSE(discretization.isTimeStageAfterImpulse(i));
    EXPECT_FALSE(discretization.isTimeStageAfterLift(i));
  }
}


TEST_P(TimeDiscretizationTest, discretizeGridBased) {
  TimeDiscretization discretization(T, N, max_num_events);
  const auto robot = GetParam();
  const auto contact_sequence = createContactSequence(robot);
  discretization.discretize(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N);
  EXPECT_EQ(discretization.N_impulse(), contact_sequence->numImpulseEvents());
  EXPECT_EQ(discretization.N_lift(), contact_sequence->numLiftEvents());
  std::vector<double> t_impulse, t_lift;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    t_impulse.push_back(contact_sequence->impulseTime(i));
    time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    t_lift.push_back(contact_sequence->liftTime(i));
    time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    EXPECT_EQ(time_stage_before_impulse[i], discretization.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], discretization.t_impulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt-t, discretization.dt(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(discretization.dt(time_stage_before_impulse[i])+discretization.dt_aux(i), dt);
    EXPECT_FALSE(discretization.isSTOEnabledImpulse(i));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    EXPECT_EQ(time_stage_before_lift[i], discretization.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(t_lift[i], discretization.t_lift(i));
    EXPECT_DOUBLE_EQ(t_lift[i]-time_stage_before_lift[i]*dt-t, discretization.dt(time_stage_before_lift[i]));
    EXPECT_DOUBLE_EQ(discretization.dt(time_stage_before_lift[i])+discretization.dt_lift(i), dt);
    EXPECT_FALSE(discretization.isSTOEnabledLift(i));
  }
  std::vector<int> time_stage_before_events;
  for (auto e :  time_stage_before_impulse) {
    time_stage_before_events.push_back(e);
  }
  for (auto e :  time_stage_before_lift) {
    time_stage_before_events.push_back(e);
  }
  std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
  time_stage_before_events.push_back(N+1);
  int contact_phase_ref = 0;
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(discretization.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<discretization.N_impulse(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeImpulse(i)+1, 
              discretization.timeStageAfterImpulse(i));
    EXPECT_EQ(discretization.contactPhaseAfterImpulse(i), 
              discretization.contactPhase(discretization.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<discretization.N_lift(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeLift(i)+1, 
              discretization.timeStageAfterLift(i));
    EXPECT_EQ(discretization.contactPhaseAfterLift(i), 
              discretization.contactPhase(discretization.timeStageAfterLift(i)));
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(discretization.timeStageBeforeImpulse(discretization.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(discretization.timeStageBeforeLift(discretization.liftIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.liftIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_DOUBLE_EQ(discretization.t(i), t+i*dt);
  }
  for (int i=0; i<N; ++i) {
    if (!discretization.isTimeStageBeforeImpulse(i) && !discretization.isTimeStageBeforeLift(i)) {
      EXPECT_DOUBLE_EQ(discretization.dt(i), dt);
    }
  }
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events; ++event_index) {
    EXPECT_FALSE(discretization.eventType(event_index)==DiscreteEventType::None);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      EXPECT_EQ(discretization.eventIndexImpulse(impulse_index), event_index);
      ++impulse_index;
    }
    else {
      EXPECT_EQ(discretization.eventIndexLift(lift_index), event_index);
      ++lift_index;
    }
  }
  EXPECT_NO_THROW(
    std::cout << discretization << std::endl;
  );
}


TEST_P(TimeDiscretizationTest, discretizeGridBased_switchingTimesOnGrids) {
  TimeDiscretization discretization(T, N, max_num_events);
  const auto robot = GetParam();
  const auto contact_sequence = createContactSequenceOnGrid(robot);
  discretization.discretize(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N-max_num_events);
  EXPECT_EQ(discretization.N_impulse(), contact_sequence->numImpulseEvents());
  EXPECT_EQ(discretization.N_lift(), contact_sequence->numLiftEvents());
  double ti = t;
  for (int i=0; i<discretization.N(); ++i) {
    EXPECT_NEAR(discretization.dt(i), dt, min_dt);
    EXPECT_NEAR(discretization.t(i), ti, min_dt);
    ti += dt;
    if (discretization.isTimeStageBeforeImpulse(i) || discretization.isTimeStageBeforeLift(i)) {
      ti += dt;
    }
  }
  EXPECT_DOUBLE_EQ(discretization.t(discretization.N()), t+T);
  EXPECT_NO_THROW(
    std::cout << discretization << std::endl;
  );
}


TEST_P(TimeDiscretizationTest, discretizeGridBased_eventTimesAreLargerThanHorizon) {
  const double T_short = 0.5 * T;
  const double dt_short = T_short / N;
  TimeDiscretization discretization(T_short, N, max_num_events);
  const auto robot = GetParam();
  const auto contact_sequence = createContactSequence(robot);
  discretization.discretize(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N);
  EXPECT_TRUE(discretization.N_impulse() <= contact_sequence->numImpulseEvents());
  EXPECT_TRUE(discretization.N_lift() <= contact_sequence->numLiftEvents());
  std::vector<double> t_impulse, t_lift;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    if (contact_sequence->impulseTime(i) > t+T_short-min_dt) {
      break;
    }
    t_impulse.push_back(contact_sequence->impulseTime(i));
    time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt_short));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    if (contact_sequence->liftTime(i) > t+T_short-min_dt) {
      break;
    }
    t_lift.push_back(contact_sequence->liftTime(i));
    time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt_short));
  }
  const int N_impulse = t_impulse.size();
  const int N_lift = t_lift.size();
  EXPECT_EQ(discretization.N_impulse(), N_impulse);
  EXPECT_EQ(discretization.N_lift(), N_lift);
  for (int i=0; i<N_impulse; ++i) {
    EXPECT_EQ(time_stage_before_impulse[i], discretization.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], discretization.t_impulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt_short-t, discretization.dt(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(discretization.dt(time_stage_before_impulse[i])+discretization.dt_aux(i), dt_short);
    EXPECT_FALSE(discretization.isSTOEnabledImpulse(i));
  }
  for (int i=0; i<N_lift; ++i) {
    EXPECT_EQ(time_stage_before_lift[i], discretization.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(t_lift[i], discretization.t_lift(i));
    EXPECT_DOUBLE_EQ(t_lift[i]-time_stage_before_lift[i]*dt_short-t, discretization.dt(time_stage_before_lift[i]));
    EXPECT_DOUBLE_EQ(discretization.dt(time_stage_before_lift[i])+discretization.dt_lift(i), dt_short);
    EXPECT_FALSE(discretization.isSTOEnabledLift(i));
  }
  std::vector<int> time_stage_before_events;
  for (auto e :  time_stage_before_impulse) {
    time_stage_before_events.push_back(e);
  }
  for (auto e :  time_stage_before_lift) {
    time_stage_before_events.push_back(e);
  }
  std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
  time_stage_before_events.push_back(N+1);
  int contact_phase_ref = 0;
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(discretization.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<discretization.N_impulse(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeImpulse(i)+1, 
              discretization.timeStageAfterImpulse(i));
    EXPECT_EQ(discretization.contactPhaseAfterImpulse(i), 
              discretization.contactPhase(discretization.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<discretization.N_lift(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeLift(i)+1, 
              discretization.timeStageAfterLift(i));
    EXPECT_EQ(discretization.contactPhaseAfterLift(i), 
              discretization.contactPhase(discretization.timeStageAfterLift(i)));
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(discretization.timeStageBeforeImpulse(discretization.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(discretization.timeStageBeforeLift(discretization.liftIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.liftIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_DOUBLE_EQ(discretization.t(i), t+i*dt_short);
  }
  for (int i=0; i<N; ++i) {
    if (!discretization.isTimeStageBeforeImpulse(i) && !discretization.isTimeStageBeforeLift(i)) {
      EXPECT_DOUBLE_EQ(discretization.dt(i), dt_short);
    }
  }
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events; ++event_index) {
    EXPECT_FALSE(discretization.eventType(event_index)==DiscreteEventType::None);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      EXPECT_EQ(discretization.eventIndexImpulse(impulse_index), event_index);
      ++impulse_index;
    }
    else {
      EXPECT_EQ(discretization.eventIndexLift(lift_index), event_index);
      ++lift_index;
    }
  }
  EXPECT_NO_THROW(
    std::cout << discretization << std::endl;
  );
  const auto time_steps = discretization.timeSteps();
  int j=0;
  for (int i=0; i<discretization.N(); ++i, ++j) {
    EXPECT_DOUBLE_EQ(time_steps[j], discretization.dt(i));
    if (discretization.isTimeStageBeforeImpulse(i)) {
      ++j;
      EXPECT_DOUBLE_EQ(time_steps[j], 
                       discretization.dt_aux(discretization.impulseIndexAfterTimeStage(i)));
    }
    else if (discretization.isTimeStageBeforeLift(i)) {
      ++j;
      EXPECT_DOUBLE_EQ(time_steps[j], 
                       discretization.dt_lift(discretization.liftIndexAfterTimeStage(i)));
    }
  }
}


TEST_P(TimeDiscretizationTest, discretizePhaseBased) {
  TimeDiscretization discretization(T, N, max_num_events);
  discretization.setDiscretizationMethod(DiscretizationMethod::PhaseBased);
  const auto robot = GetParam();
  const auto contact_sequence = createContactSequence(robot);
  discretization.meshRefinement(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N);
  EXPECT_EQ(discretization.N_impulse(), contact_sequence->numImpulseEvents());
  EXPECT_EQ(discretization.N_lift(), contact_sequence->numLiftEvents());
  std::vector<double> t_impulse, t_lift;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    t_impulse.push_back(contact_sequence->impulseTime(i));
    time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    t_lift.push_back(contact_sequence->liftTime(i));
    time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    EXPECT_EQ(time_stage_before_impulse[i], discretization.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], discretization.t_impulse(i));
    EXPECT_FALSE(discretization.isSTOEnabledImpulse(i));
    EXPECT_EQ(discretization.isSTOEnabledImpulse(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexImpulse(i)));
    EXPECT_EQ(discretization.isSTOEnabledImpulse(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexImpulse(i)));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    EXPECT_EQ(time_stage_before_lift[i], discretization.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(t_lift[i], discretization.t_lift(i));
    EXPECT_FALSE(discretization.isSTOEnabledLift(i));
    EXPECT_EQ(discretization.isSTOEnabledLift(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexLift(i)));
  }
  EXPECT_EQ(discretization.isSTOEnabledPhase(0), 
            discretization.isSTOEnabledEvent(0));
  for (int i=1; i<contact_sequence->numDiscreteEvents(); ++i) {
    EXPECT_EQ(discretization.isSTOEnabledPhase(i), 
              (discretization.isSTOEnabledEvent(i-1)||discretization.isSTOEnabledEvent(i)));
  }
  EXPECT_EQ(discretization.isSTOEnabledPhase(contact_sequence->numDiscreteEvents()), 
            discretization.isSTOEnabledEvent(contact_sequence->numDiscreteEvents()-1));
  std::vector<int> time_stage_before_events;
  for (auto e :  time_stage_before_impulse) {
    time_stage_before_events.push_back(e);
  }
  for (auto e :  time_stage_before_lift) {
    time_stage_before_events.push_back(e);
  }
  std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
  time_stage_before_events.push_back(N+1);
  int contact_phase_ref = 0;
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(discretization.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<discretization.N_impulse(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeImpulse(i)+1, 
              discretization.timeStageAfterImpulse(i));
    EXPECT_EQ(discretization.contactPhaseAfterImpulse(i), 
              discretization.contactPhase(discretization.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<discretization.N_lift(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeLift(i)+1, 
              discretization.timeStageAfterLift(i));
    EXPECT_EQ(discretization.contactPhaseAfterLift(i), 
              discretization.contactPhase(discretization.timeStageAfterLift(i)));
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(discretization.timeStageBeforeImpulse(discretization.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(discretization.timeStageBeforeLift(discretization.liftIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.liftIndexAfterTimeStage(i), -1);
    }
  }
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  int time_stage_before_event = 0;
  double t_prev_event = t;
  double dt_prev_aux = discretization.dt(0);
  for (int event_index=0; event_index<num_events; ++event_index) {
    ASSERT_FALSE(discretization.eventType(event_index)==DiscreteEventType::None);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      EXPECT_EQ(discretization.eventIndexImpulse(impulse_index), event_index);
      const int grids_phase = discretization.timeStageBeforeImpulse(impulse_index) 
                              - time_stage_before_event + 1;
      EXPECT_EQ(discretization.N_phase(event_index), grids_phase);
      const double dt_phase = (discretization.t_impulse(impulse_index)-t_prev_event) / grids_phase;
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeImpulse(impulse_index); ++stage) {
        EXPECT_DOUBLE_EQ(discretization.dt(stage), dt_phase);
      }
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeImpulse(impulse_index)-1; ++stage) {
        EXPECT_NEAR(discretization.t(stage+1)-discretization.t(stage), dt_phase, min_dt);
      }
      EXPECT_NEAR(discretization.t_impulse(impulse_index)-discretization.t(discretization.timeStageBeforeImpulse(impulse_index)), 
                  dt_phase, min_dt);
      EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
      time_stage_before_event = discretization.timeStageBeforeImpulse(impulse_index);
      t_prev_event = discretization.t_impulse(impulse_index);
      dt_prev_aux = discretization.dt_aux(impulse_index);
      ++impulse_index;
    }
    else {
      EXPECT_EQ(discretization.eventIndexLift(lift_index), event_index);
      const int grids_phase = discretization.timeStageBeforeLift(lift_index) 
                              - time_stage_before_event + 1;
      EXPECT_EQ(discretization.N_phase(event_index), grids_phase);
      const double dt_phase = (discretization.t_lift(lift_index)-t_prev_event) / grids_phase;
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeLift(lift_index); ++stage) {
        EXPECT_DOUBLE_EQ(discretization.dt(stage), dt_phase);
      }
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeLift(lift_index)-1; ++stage) {
        EXPECT_NEAR(discretization.t(stage+1)-discretization.t(stage), dt_phase, min_dt);
      }
      EXPECT_NEAR(discretization.t_lift(lift_index)-discretization.t(discretization.timeStageBeforeLift(lift_index)), 
                  dt_phase, min_dt);
      EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
      time_stage_before_event = discretization.timeStageBeforeLift(lift_index);
      t_prev_event = discretization.t_lift(lift_index);
      dt_prev_aux = discretization.dt_lift(lift_index);
      ++lift_index;
    }
  }
  EXPECT_NO_THROW(
    std::cout << discretization << std::endl;
  );
}


TEST_P(TimeDiscretizationTest, discretizePhaseBased_eventTimesAreLargerThanHorizon) {
  const double T_short = 0.5 * T;
  const double dt_short = T_short / N;
  TimeDiscretization discretization(T_short, N, max_num_events);
  discretization.setDiscretizationMethod(DiscretizationMethod::PhaseBased);
  const auto robot = GetParam();
  const auto contact_sequence = createContactSequence(robot);
  discretization.meshRefinement(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N);
  std::vector<double> t_impulse, t_lift;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
    if (contact_sequence->impulseTime(i) >= t+T_short-min_dt) {
      break;
    }
    t_impulse.push_back(contact_sequence->impulseTime(i));
    time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt_short));
  }
  for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
    if (contact_sequence->liftTime(i) >= t+T_short-min_dt) {
      break;
    }
    t_lift.push_back(contact_sequence->liftTime(i));
    time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt_short));
  }
  const int N_impulse = t_impulse.size();
  const int N_lift = t_lift.size();
  EXPECT_EQ(discretization.N_impulse(), N_impulse);
  EXPECT_EQ(discretization.N_lift(), N_lift);
  for (int i=0; i<N_impulse; ++i) {
    EXPECT_EQ(time_stage_before_impulse[i], discretization.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], discretization.t_impulse(i));
    EXPECT_FALSE(discretization.isSTOEnabledImpulse(i));
    EXPECT_EQ(discretization.isSTOEnabledImpulse(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexImpulse(i)));
    EXPECT_EQ(discretization.isSTOEnabledImpulse(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexImpulse(i)));
  }
  for (int i=0; i<N_lift; ++i) {
    EXPECT_EQ(time_stage_before_lift[i], discretization.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(t_lift[i], discretization.t_lift(i));
    EXPECT_FALSE(discretization.isSTOEnabledLift(i));
    EXPECT_EQ(discretization.isSTOEnabledLift(i), 
              discretization.isSTOEnabledEvent(discretization.eventIndexLift(i)));
  }
  EXPECT_EQ(discretization.isSTOEnabledPhase(0), 
            discretization.isSTOEnabledEvent(0));
  for (int i=1; i<N_lift+N_impulse; ++i) {
    EXPECT_EQ(discretization.isSTOEnabledPhase(i), 
              (discretization.isSTOEnabledEvent(i-1)||discretization.isSTOEnabledEvent(i)));
  }
  EXPECT_EQ(discretization.isSTOEnabledPhase(N_lift+N_impulse), 
            discretization.isSTOEnabledEvent(N_lift+N_impulse-1));
  std::vector<int> time_stage_before_events;
  for (auto e :  time_stage_before_impulse) {
    time_stage_before_events.push_back(e);
  }
  for (auto e :  time_stage_before_lift) {
    time_stage_before_events.push_back(e);
  }
  std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
  time_stage_before_events.push_back(N+1);
  int contact_phase_ref = 0;
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(discretization.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<discretization.N_impulse(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeImpulse(i)+1, 
              discretization.timeStageAfterImpulse(i));
    EXPECT_EQ(discretization.contactPhaseAfterImpulse(i), 
              discretization.contactPhase(discretization.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<discretization.N_lift(); ++i) {
    EXPECT_EQ(discretization.timeStageBeforeLift(i)+1, 
              discretization.timeStageAfterLift(i));
    EXPECT_EQ(discretization.contactPhaseAfterLift(i), 
              discretization.contactPhase(discretization.timeStageAfterLift(i)));
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(discretization.timeStageBeforeImpulse(discretization.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (discretization.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(discretization.timeStageBeforeLift(discretization.liftIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(discretization.liftIndexAfterTimeStage(i), -1);
    }
  }
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  int time_stage_before_event = 0;
  double t_prev_event = t;
  double dt_prev_aux = discretization.dt(0);
  for (int event_index=0; event_index<num_events; ++event_index) {
    ASSERT_FALSE(discretization.eventType(event_index)==DiscreteEventType::None);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      EXPECT_EQ(discretization.eventIndexImpulse(impulse_index), event_index);
      const int grids_phase = discretization.timeStageBeforeImpulse(impulse_index) 
                              - time_stage_before_event + 1;
      EXPECT_EQ(discretization.N_phase(event_index), grids_phase);
      const double dt_phase = (discretization.t_impulse(impulse_index)-t_prev_event) / grids_phase;
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeImpulse(impulse_index); ++stage) {
        EXPECT_DOUBLE_EQ(discretization.dt(stage), dt_phase);
      }
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeImpulse(impulse_index)-1; ++stage) {
        EXPECT_NEAR(discretization.t(stage+1)-discretization.t(stage), dt_phase, min_dt);
      }
      EXPECT_NEAR(discretization.t_impulse(impulse_index)-discretization.t(discretization.timeStageBeforeImpulse(impulse_index)), 
                  dt_phase, min_dt);
      EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
      time_stage_before_event = discretization.timeStageBeforeImpulse(impulse_index);
      t_prev_event = discretization.t_impulse(impulse_index);
      dt_prev_aux = discretization.dt_aux(impulse_index);
      ++impulse_index;
    }
    else {
      EXPECT_EQ(discretization.eventIndexLift(lift_index), event_index);
      const int grids_phase = discretization.timeStageBeforeLift(lift_index) 
                              - time_stage_before_event + 1;
      EXPECT_EQ(discretization.N_phase(event_index), grids_phase);
      const double dt_phase = (discretization.t_lift(lift_index)-t_prev_event) / grids_phase;
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeLift(lift_index); ++stage) {
        EXPECT_DOUBLE_EQ(discretization.dt(stage), dt_phase);
      }
      for (int stage=time_stage_before_event+1; 
            stage<discretization.timeStageBeforeLift(lift_index)-1; ++stage) {
        EXPECT_NEAR(discretization.t(stage+1)-discretization.t(stage), dt_phase, min_dt);
      }
      EXPECT_NEAR(discretization.t_lift(lift_index)-discretization.t(discretization.timeStageBeforeLift(lift_index)), 
                  dt_phase, min_dt);
      EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
      time_stage_before_event = discretization.timeStageBeforeLift(lift_index);
      t_prev_event = discretization.t_lift(lift_index);
      dt_prev_aux = discretization.dt_lift(lift_index);
      ++lift_index;
    }
  }
  EXPECT_NO_THROW(
    std::cout << discretization << std::endl;
  );
  const auto time_steps = discretization.timeSteps();
  int j=0;
  for (int i=0; i<discretization.N(); ++i, ++j) {
    EXPECT_DOUBLE_EQ(time_steps[j], discretization.dt(i));
    if (discretization.isTimeStageBeforeImpulse(i)) {
      ++j;
      EXPECT_DOUBLE_EQ(time_steps[j], 
                       discretization.dt_aux(discretization.impulseIndexAfterTimeStage(i)));
    }
    else if (discretization.isTimeStageBeforeLift(i)) {
      ++j;
      EXPECT_DOUBLE_EQ(time_steps[j], 
                       discretization.dt_lift(discretization.liftIndexAfterTimeStage(i)));
    }
  }
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, TimeDiscretizationTest, 
  ::testing::Values(testhelper::CreateFixedBaseRobot(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateFloatingBaseRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}