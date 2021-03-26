#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/robot/robot.hpp"

#include "robot_factory.hpp"


namespace idocp {

class OCPDiscretizerTest : public ::testing::Test {
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

  void testConstructor(const Robot& robot) const;
  void testDiscretizeOCP(const Robot& robot) const;
  void testDiscretizeOCPOnGrid(const Robot& robot) const;

  int N, max_num_events;
  double t, T, dt, min_dt;
};


ContactSequence OCPDiscretizerTest::createContactSequence(const Robot& robot) const {
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
    contact_sequence.push_back(tmp, event_time);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


ContactSequence OCPDiscretizerTest::createContactSequenceOnGrid(const Robot& robot) const {
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
    contact_sequence.push_back(tmp, event_time);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


void OCPDiscretizerTest::testConstructor(const Robot& robot) const {
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  EXPECT_EQ(ocp_discretizer.N(), N);
  EXPECT_EQ(ocp_discretizer.N_impulse(), 0);
  EXPECT_EQ(ocp_discretizer.N_lift(), 0);
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(ocp_discretizer.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(ocp_discretizer.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(ocp_discretizer.isTimeStageBeforeLift(i));
  }
  for (int i=1; i<=N; ++i) {
    EXPECT_FALSE(ocp_discretizer.isTimeStageAfterImpulse(i));
    EXPECT_FALSE(ocp_discretizer.isTimeStageAfterLift(i));
  }
}


void OCPDiscretizerTest::testDiscretizeOCP(const Robot& robot) const {
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(ocp_discretizer.N(), N);
  EXPECT_EQ(ocp_discretizer.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(ocp_discretizer.N_lift(), contact_sequence.numLiftEvents());
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
    EXPECT_EQ(time_stage_before_impulse[i], ocp_discretizer.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i], ocp_discretizer.t_impulse(i));
    EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt-t, ocp_discretizer.dt(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(ocp_discretizer.dt(time_stage_before_impulse[i])+ocp_discretizer.dt_aux(i), dt);
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    EXPECT_EQ(time_stage_before_lift[i], ocp_discretizer.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(t_lift[i], ocp_discretizer.t_lift(i));
    EXPECT_DOUBLE_EQ(t_lift[i]-time_stage_before_lift[i]*dt-t, ocp_discretizer.dt(time_stage_before_lift[i]));
    EXPECT_DOUBLE_EQ(ocp_discretizer.dt(time_stage_before_lift[i])+ocp_discretizer.dt_lift(i), dt);
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
    EXPECT_EQ(ocp_discretizer.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<ocp_discretizer.N_impulse(); ++i) {
    EXPECT_EQ(ocp_discretizer.timeStageBeforeImpulse(i)+1, 
              ocp_discretizer.timeStageAfterImpulse(i));
    EXPECT_EQ(ocp_discretizer.contactPhaseAfterImpulse(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<ocp_discretizer.N_lift(); ++i) {
    EXPECT_EQ(ocp_discretizer.timeStageBeforeLift(i)+1, 
              ocp_discretizer.timeStageAfterLift(i));
    EXPECT_EQ(ocp_discretizer.contactPhaseAfterLift(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageAfterLift(i)));
  }
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(ocp_discretizer.timeStageBeforeImpulse(ocp_discretizer.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(ocp_discretizer.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(ocp_discretizer.timeStageBeforeLift(ocp_discretizer.liftIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(ocp_discretizer.liftIndexAfterTimeStage(i), -1);
    }
  }
  const double dt = T/N;
  for (int i=0; i<=N; ++i) {
    EXPECT_DOUBLE_EQ(ocp_discretizer.t(i), t+i*dt);
  }
  for (int i=0; i<N; ++i) {
    if (!ocp_discretizer.isTimeStageBeforeImpulse(i) && !ocp_discretizer.isTimeStageBeforeLift(i)) {
      EXPECT_DOUBLE_EQ(ocp_discretizer.dt(i), dt);
    }
  }
}


void OCPDiscretizerTest::testDiscretizeOCPOnGrid(const Robot& robot) const {
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequenceOnGrid(robot);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(ocp_discretizer.N(), N-max_num_events);
  EXPECT_EQ(ocp_discretizer.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(ocp_discretizer.N_lift(), contact_sequence.numLiftEvents());
  double ti = t;
  for (int i=0; i<ocp_discretizer.N(); ++i) {
    EXPECT_NEAR(ocp_discretizer.dt(i), dt, min_dt);
    EXPECT_NEAR(ocp_discretizer.t(i), ti, min_dt);
    ti += dt;
    if (ocp_discretizer.isTimeStageBeforeImpulse(i) || ocp_discretizer.isTimeStageBeforeLift(i)) {
      ti += dt;
    }
  }
  EXPECT_DOUBLE_EQ(ocp_discretizer.t(ocp_discretizer.N()), t+T);
}


TEST_F(OCPDiscretizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testConstructor(robot);
  testDiscretizeOCP(robot);
  testDiscretizeOCPOnGrid(robot);
}


TEST_F(OCPDiscretizerTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testConstructor(robot);
  testDiscretizeOCP(robot);
  testDiscretizeOCPOnGrid(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}