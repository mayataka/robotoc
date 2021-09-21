#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_ocp_discretization.hpp"
#include "idocp/robot/robot.hpp"

#include "robot_factory.hpp"


namespace idocp {

class HybridOCPDiscretizationTest : public ::testing::Test {
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

  ContactSequence createContactSequence(const Robot& robot) const;
  ContactSequence createContactSequenceOnGrid(const Robot& robot) const;

  void test_constructor(const Robot& robot) const;
  void test_discretizeOCP(const Robot& robot) const;
  void test_discretizeOCPOnGrid(const Robot& robot) const;

  int N, max_num_events;
  double t, T, dt, min_dt;
};


ContactSequence HybridOCPDiscretizationTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, max_num_events);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot.maxPointContacts());
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_time = t + i * event_period + dt * std::abs(Eigen::VectorXd::Random(1)[0]);
    contact_sequence.push_back(tmp, event_time, false);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


ContactSequence HybridOCPDiscretizationTest::createContactSequenceOnGrid(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, max_num_events);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot.maxPointContacts());
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_time = t + (i+1) * event_period + min_dt * Eigen::VectorXd::Random(1)[0];
    contact_sequence.push_back(tmp, event_time, false);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


void HybridOCPDiscretizationTest::test_constructor(const Robot& robot) const {
  HybridOCPDiscretization discretization(T, N, max_num_events);
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


void HybridOCPDiscretizationTest::test_discretizeOCP(const Robot& robot) const {
  HybridOCPDiscretization discretization(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  discretization.discretize(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N);
  EXPECT_EQ(discretization.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(discretization.N_lift(), contact_sequence.numLiftEvents());
  std::vector<double> t_impulse, t_lift;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    t_impulse.push_back(contact_sequence.impulseTime(i));
    time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    t_lift.push_back(contact_sequence.liftTime(i));
    time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt));
  }
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    EXPECT_EQ(time_stage_before_impulse[i], discretization.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], discretization.t_impulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt-t, discretization.dt(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(discretization.dt(time_stage_before_impulse[i])+discretization.dt_aux(i), dt);
    EXPECT_FALSE(discretization.isSTOEnabledImpulse(i));
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
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
  const double dt = T/N;
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
  EXPECT_NO_THROW(discretization.showInfo());
}


void HybridOCPDiscretizationTest::test_discretizeOCPOnGrid(const Robot& robot) const {
  HybridOCPDiscretization discretization(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequenceOnGrid(robot);
  discretization.discretize(contact_sequence, t);
  EXPECT_EQ(discretization.N(), N-max_num_events);
  EXPECT_EQ(discretization.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(discretization.N_lift(), contact_sequence.numLiftEvents());
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
  EXPECT_NO_THROW(discretization.showInfo());
}


TEST_F(HybridOCPDiscretizationTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test_constructor(robot);
  test_discretizeOCP(robot);
  test_discretizeOCPOnGrid(robot);
}


TEST_F(HybridOCPDiscretizationTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test_constructor(robot);
  test_discretizeOCP(robot);
  test_discretizeOCPOnGrid(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}