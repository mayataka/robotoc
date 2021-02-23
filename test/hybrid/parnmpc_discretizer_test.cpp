#include <random>
#include <vector>
#include <algorithm>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"
#include "idocp/robot/robot.hpp"

namespace idocp {

class ParNMPCDiscretizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
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

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_events;
  double t, T, dt, min_dt;
};


ContactSequence ParNMPCDiscretizerTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, N);
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


ContactSequence ParNMPCDiscretizerTest::createContactSequenceOnGrid(const Robot& robot) const {
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


void ParNMPCDiscretizerTest::testConstructor(const Robot& robot) const {
  ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_events);
  EXPECT_EQ(parnmpc_discretizer.N(), N);
  EXPECT_EQ(parnmpc_discretizer.N_impulse(), 0);
  EXPECT_EQ(parnmpc_discretizer.N_lift(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(parnmpc_discretizer.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageAfterImpulse(i));
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageAfterLift(i));
  }
  for (int i=0; i<N-1; ++i) {
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageBeforeLift(i));
  }
}

void ParNMPCDiscretizerTest::testDiscretizeOCP(const Robot& robot) const {
  ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  parnmpc_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(parnmpc_discretizer.N(), N);
  EXPECT_EQ(parnmpc_discretizer.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(parnmpc_discretizer.N_lift(), contact_sequence.numLiftEvents());
  std::vector<double> impulse_time, lift_time;
  std::vector<int> time_stage_after_impulse, time_stage_after_lift;
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    impulse_time.push_back(contact_sequence.impulseTime(i));
    time_stage_after_impulse.push_back(std::ceil((impulse_time[i]-t)/dt)-1);
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    lift_time.push_back(contact_sequence.liftTime(i));
    time_stage_after_lift.push_back(std::ceil((lift_time[i]-t)/dt)-1);
  }
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    EXPECT_EQ(time_stage_after_impulse[i], parnmpc_discretizer.timeStageAfterImpulse(i));
    EXPECT_EQ(parnmpc_discretizer.timeStageAfterImpulse(i), parnmpc_discretizer.timeStageBeforeImpulse(i)+1);
    EXPECT_DOUBLE_EQ(impulse_time[i], parnmpc_discretizer.t_impulse(i));
    EXPECT_DOUBLE_EQ((time_stage_after_impulse[i]+1)*dt+t-impulse_time[i], parnmpc_discretizer.dt(time_stage_after_impulse[i]));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.dt(time_stage_after_impulse[i])+parnmpc_discretizer.dt_aux(i), dt);
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    EXPECT_EQ(time_stage_after_lift[i], parnmpc_discretizer.timeStageAfterLift(i));
    EXPECT_EQ(parnmpc_discretizer.timeStageAfterLift(i), parnmpc_discretizer.timeStageBeforeLift(i)+1);
    EXPECT_DOUBLE_EQ(lift_time[i], parnmpc_discretizer.t_lift(i));
    EXPECT_DOUBLE_EQ((time_stage_after_lift[i]+1)*dt+t-lift_time[i], parnmpc_discretizer.dt(time_stage_after_lift[i]));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.dt(time_stage_after_lift[i])+parnmpc_discretizer.dt_lift(i), dt);
  }
  std::vector<int> time_stage_after_events;
  for (auto e :  time_stage_after_impulse) {
    time_stage_after_events.push_back(e);
  }
  for (auto e :  time_stage_after_lift) {
    time_stage_after_events.push_back(e);
  }
  std::sort(time_stage_after_events.begin(), time_stage_after_events.end());
  time_stage_after_events.push_back(N);
  int contact_phase_ref = 0;
  for (int i=0; i<N; ++i) {
    if (i == time_stage_after_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
    EXPECT_EQ(parnmpc_discretizer.contactPhase(i), contact_phase_ref);
  }
  for (int i=0; i<parnmpc_discretizer.N_impulse(); ++i) {
    EXPECT_EQ(parnmpc_discretizer.timeStageBeforeImpulse(i)+1, 
              parnmpc_discretizer.timeStageAfterImpulse(i));
    if (parnmpc_discretizer.timeStageBeforeImpulse(i) >= 0) {
      EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeImpulse(i), 
                parnmpc_discretizer.contactPhase(parnmpc_discretizer.timeStageBeforeImpulse(i)));
    }
    else {
      EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeImpulse(i), 0);
    }
  }
  for (int i=0; i<parnmpc_discretizer.N_lift(); ++i) {
    EXPECT_EQ(parnmpc_discretizer.timeStageBeforeLift(i)+1, 
              parnmpc_discretizer.timeStageAfterLift(i));
    if (parnmpc_discretizer.timeStageBeforeLift(i) >= 0) {
      EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeLift(i), 
                parnmpc_discretizer.contactPhase(parnmpc_discretizer.timeStageBeforeLift(i)));
    }
    else {
      EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeLift(i), 0);
    }
  }
  for (int i=0; i<N; ++i) {
    if (parnmpc_discretizer.isTimeStageAfterImpulse(i)) {
      EXPECT_EQ(parnmpc_discretizer.timeStageAfterImpulse(parnmpc_discretizer.impulseIndexBeforeTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(parnmpc_discretizer.impulseIndexBeforeTimeStage(i), -1);
    }
  }
  for (int i=0; i<N; ++i) {
    if (parnmpc_discretizer.isTimeStageAfterLift(i)) {
      EXPECT_EQ(parnmpc_discretizer.timeStageAfterLift(parnmpc_discretizer.liftIndexBeforeTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(parnmpc_discretizer.liftIndexBeforeTimeStage(i), -1);
    }
  }
  const double dt = T/N;
  for (int i=0; i<N; ++i) {
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.t(i), t+(i+1)*dt);
  }
  for (int i=0; i<N; ++i) {
    if (!parnmpc_discretizer.isTimeStageAfterImpulse(i) && !parnmpc_discretizer.isTimeStageAfterLift(i)) {
      EXPECT_DOUBLE_EQ(parnmpc_discretizer.dt(i), dt);
    }
  }
}


void ParNMPCDiscretizerTest::testDiscretizeOCPOnGrid(const Robot& robot) const {
  ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequenceOnGrid(robot);
  parnmpc_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(parnmpc_discretizer.N(), N-max_num_events);
  EXPECT_EQ(parnmpc_discretizer.N_impulse(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(parnmpc_discretizer.N_lift(), contact_sequence.numLiftEvents());
  double ti = t+dt;
  for (int i=0; i<parnmpc_discretizer.N(); ++i) {
    EXPECT_NEAR(parnmpc_discretizer.dt(i), dt, min_dt);
    EXPECT_NEAR(parnmpc_discretizer.t(i), ti, min_dt);
    ti += dt;
    if (parnmpc_discretizer.isTimeStageBeforeImpulse(i) || parnmpc_discretizer.isTimeStageBeforeLift(i)) {
      ti += dt;
    }
  }
}


TEST_F(ParNMPCDiscretizerTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testDiscretizeOCP(robot);
}


TEST_F(ParNMPCDiscretizerTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testDiscretizeOCP(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}