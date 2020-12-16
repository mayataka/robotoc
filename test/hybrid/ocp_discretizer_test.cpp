#include <random>
#include <vector>
#include <algorithm>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/robot/robot.hpp"

namespace idocp {

class OCPDiscretizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_events = 5;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    T = 1;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;

  void testConstructor(const Robot& robot) const;
  void testDiscretizeOCP(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_events;
  double t, T, dtau;
};


ContactSequence OCPDiscretizerTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, N);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_period = 3 * dtau;
    tmp.eventTime = t + i * event_period + 0.1 * event_period * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_events; ++i) {
    contact_sequence.pushBackDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void OCPDiscretizerTest::testConstructor(const Robot& robot) const {
  OCPDiscretizer ocp_discrertizer(T, N, max_num_events);
  EXPECT_EQ(ocp_discrertizer.N(), N);
  EXPECT_EQ(ocp_discrertizer.numImpulseStages(), 0);
  EXPECT_EQ(ocp_discrertizer.numLiftStages(), 0);
  EXPECT_FALSE(ocp_discrertizer.existImpulse());
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(ocp_discrertizer.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(ocp_discrertizer.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(ocp_discrertizer.isTimeStageBeforeLift(i));
  }
}

void OCPDiscretizerTest::testDiscretizeOCP(const Robot& robot) const {
  OCPDiscretizer ocp_discrertizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  ocp_discrertizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(ocp_discrertizer.N(), N);
  EXPECT_EQ(ocp_discrertizer.numImpulseStages(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(ocp_discrertizer.numLiftStages(), contact_sequence.numLiftEvents());
  EXPECT_EQ(ocp_discrertizer.existImpulse(), (contact_sequence.numImpulseEvents() > 0));
  std::vector<double> impulse_time, lift_time;
  std::vector<int> time_stage_before_impulse, time_stage_before_lift;
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    impulse_time.push_back(contact_sequence.impulseTime(i)-t);
    time_stage_before_impulse.push_back(std::floor(impulse_time[i]/dtau));
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    lift_time.push_back(contact_sequence.liftTime(i)-t);
    time_stage_before_lift.push_back(std::floor(lift_time[i]/dtau));
  }
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    EXPECT_DOUBLE_EQ(time_stage_before_impulse[i], ocp_discrertizer.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(impulse_time[i], ocp_discrertizer.t_impulse(i));
    EXPECT_DOUBLE_EQ(impulse_time[i]-time_stage_before_impulse[i]*dtau, ocp_discrertizer.dtau(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(ocp_discrertizer.dtau(time_stage_before_impulse[i])+ocp_discrertizer.dtau_aux(i), dtau);
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    EXPECT_DOUBLE_EQ(time_stage_before_lift[i], ocp_discrertizer.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(lift_time[i], ocp_discrertizer.t_lift(i));
    EXPECT_DOUBLE_EQ(lift_time[i]-time_stage_before_lift[i]*dtau, ocp_discrertizer.dtau(time_stage_before_lift[i]));
    EXPECT_DOUBLE_EQ(ocp_discrertizer.dtau(time_stage_before_lift[i])+ocp_discrertizer.dtau_lift(i), dtau);
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
    EXPECT_EQ(ocp_discrertizer.contactPhase(i), contact_phase_ref);
    if (i == time_stage_before_events[contact_phase_ref]) {
      ++contact_phase_ref;
    }
  }
  for (int i=0; i<ocp_discrertizer.numImpulseStages(); ++i) {
    EXPECT_EQ(ocp_discrertizer.timeStageBeforeImpulse(i)+1, 
              ocp_discrertizer.timeStageAfterImpulse(i));
  }
  for (int i=0; i<ocp_discrertizer.numLiftStages(); ++i) {
    EXPECT_EQ(ocp_discrertizer.timeStageBeforeLift(i)+1, 
              ocp_discrertizer.timeStageAfterLift(i));
  }
  for (int i=0; i<=N; ++i) {
    if (ocp_discrertizer.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(ocp_discrertizer.timeStageBeforeImpulse(ocp_discrertizer.impulseIndex(i)), i);
    }
    else {
      EXPECT_EQ(ocp_discrertizer.impulseIndex(i), -1);
    }
  }
  for (int i=0; i<=N; ++i) {
    if (ocp_discrertizer.isTimeStageBeforeLift(i)) {
      EXPECT_EQ(ocp_discrertizer.timeStageBeforeLift(ocp_discrertizer.liftIndex(i)), i);
    }
    else {
      EXPECT_EQ(ocp_discrertizer.liftIndex(i), -1);
    }
  }
}


TEST_F(OCPDiscretizerTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testDiscretizeOCP(robot);
}


TEST_F(OCPDiscretizerTest, floatingBase) {
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