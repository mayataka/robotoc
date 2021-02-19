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
    DiscreteEvent tmp(robot.maxPointContacts());
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_period = 3 * dtau;
    const double event_time = t + i * event_period + 0.1 * event_period * std::abs(Eigen::VectorXd::Random(1)[0]);
    contact_sequence.push_back(tmp, event_time);
    pre_contact_status = post_contact_status;
  }
  return contact_sequence;
}


void OCPDiscretizerTest::testConstructor(const Robot& robot) const {
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  EXPECT_EQ(ocp_discretizer.N(), N);
  EXPECT_EQ(ocp_discretizer.numImpulseStages(), 0);
  EXPECT_EQ(ocp_discretizer.numLiftStages(), 0);
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(ocp_discretizer.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(ocp_discretizer.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(ocp_discretizer.isTimeStageBeforeLift(i));
  }
}

void OCPDiscretizerTest::testDiscretizeOCP(const Robot& robot) const {
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(ocp_discretizer.N(), N);
  EXPECT_EQ(ocp_discretizer.numImpulseStages(), contact_sequence.numImpulseEvents());
  EXPECT_EQ(ocp_discretizer.numLiftStages(), contact_sequence.numLiftEvents());
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
    EXPECT_EQ(time_stage_before_impulse[i], ocp_discretizer.timeStageBeforeImpulse(i));
    EXPECT_DOUBLE_EQ(impulse_time[i], ocp_discretizer.t_impulse(i));
    EXPECT_DOUBLE_EQ(impulse_time[i]-time_stage_before_impulse[i]*dtau, ocp_discretizer.dtau(time_stage_before_impulse[i]));
    EXPECT_DOUBLE_EQ(ocp_discretizer.dtau(time_stage_before_impulse[i])+ocp_discretizer.dtau_aux(i), dtau);
    // std::cout << ocp_discretizer.dtau(time_stage_before_impulse[i])+ocp_discretizer.dtau_aux(i) - dtau << std::endl;
    // std::cout << (time_stage_before_impulse[i]+1)*dtau-ocp_discretizer.t_impulse(i) - ocp_discretizer.dtau_aux(i) << std::endl;
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    EXPECT_EQ(time_stage_before_lift[i], ocp_discretizer.timeStageBeforeLift(i));
    EXPECT_DOUBLE_EQ(lift_time[i], ocp_discretizer.t_lift(i));
    EXPECT_DOUBLE_EQ(lift_time[i]-time_stage_before_lift[i]*dtau, ocp_discretizer.dtau(time_stage_before_lift[i]));
    EXPECT_DOUBLE_EQ(ocp_discretizer.dtau(time_stage_before_lift[i])+ocp_discretizer.dtau_lift(i), dtau);
    // std::cout << ocp_discretizer.dtau(time_stage_before_lift[i])+ocp_discretizer.dtau_lift(i) - dtau << std::endl;
    // std::cout << (time_stage_before_lift[i]+1)*dtau-ocp_discretizer.t_lift(i) - ocp_discretizer.dtau_lift(i) << std::endl;
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
  for (int i=0; i<ocp_discretizer.numImpulseStages(); ++i) {
    EXPECT_EQ(ocp_discretizer.timeStageBeforeImpulse(i)+1, 
              ocp_discretizer.timeStageAfterImpulse(i));
    EXPECT_EQ(ocp_discretizer.contactPhaseBeforeImpulse(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageBeforeImpulse(i)));
    EXPECT_EQ(ocp_discretizer.contactPhaseAfterImpulse(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageAfterImpulse(i)));
  }
  for (int i=0; i<ocp_discretizer.numLiftStages(); ++i) {
    EXPECT_EQ(ocp_discretizer.timeStageBeforeLift(i)+1, 
              ocp_discretizer.timeStageAfterLift(i));
    EXPECT_EQ(ocp_discretizer.contactPhaseBeforeLift(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageBeforeLift(i)));
    EXPECT_EQ(ocp_discretizer.contactPhaseAfterLift(i), 
              ocp_discretizer.contactPhase(ocp_discretizer.timeStageAfterLift(i)));
  }
  for (int i=0; i<=N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      EXPECT_EQ(ocp_discretizer.timeStageBeforeImpulse(ocp_discretizer.impulseIndexAfterTimeStage(i)), i);
    }
    else {
      EXPECT_EQ(ocp_discretizer.impulseIndexAfterTimeStage(i), -1);
    }
  }
  for (int i=0; i<=N; ++i) {
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
  for (int i=0; i<=N; ++i) {
    if (!ocp_discretizer.isTimeStageBeforeImpulse(i) && !ocp_discretizer.isTimeStageBeforeLift(i)) {
      EXPECT_DOUBLE_EQ(ocp_discretizer.dtau(i), dt);
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