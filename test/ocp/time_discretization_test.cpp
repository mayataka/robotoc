#include <random>
#include <vector>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/time_discretization.hpp"
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
  contact_sequence->init(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(pre_contact_status, post_contact_status);
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
  contact_sequence->init(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  const double event_period = 3 * dt;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(pre_contact_status, post_contact_status);
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
  TimeDiscretization time_discretization(T, N, max_num_events);
  EXPECT_EQ(time_discretization.N(), N);
  EXPECT_EQ(time_discretization.size(), N+1);
  for (int i=0; i<=N; ++i) {
    EXPECT_EQ(time_discretization[i].phase, 0);
  }
}


// TEST_P(TimeDiscretizationTest, discretizeGridBased) {
//   TimeDiscretization time_discretization(T, N, max_num_events);
//   const auto robot = GetParam();
//   const auto contact_sequence = createContactSequence(robot);
//   time_discretization.discretize(contact_sequence, t);
//   EXPECT_EQ(time_discretization.N(), N);
//   EXPECT_EQ(time_discretization.N_impulse(), contact_sequence->numImpulseEvents());
//   EXPECT_EQ(time_discretization.N_lift(), contact_sequence->numLiftEvents());
//   std::vector<double> t_impulse, t_lift;
//   std::vector<int> time_stage_before_impulse, time_stage_before_lift;
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     t_impulse.push_back(contact_sequence->impulseTime(i));
//     time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     t_lift.push_back(contact_sequence->liftTime(i));
//     time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt));
//   }
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     EXPECT_EQ(time_stage_before_impulse[i], time_discretization.timeStageBeforeImpulse(i));
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.impulseTime(i));
//     EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt-t, time_discretization.gridInfo(time_stage_before_impulse[i]).dt);
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(time_stage_before_impulse[i]).dt+time_discretization.gridInfoAux(i).dt, dt);
//     EXPECT_FALSE(time_discretization.isSTOEnabledImpulse(i));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     EXPECT_EQ(time_stage_before_lift[i], time_discretization.timeStageBeforeLift(i));
//     EXPECT_DOUBLE_EQ(t_lift[i], time_discretization.liftTime(i));
//     EXPECT_DOUBLE_EQ(t_lift[i]-time_stage_before_lift[i]*dt-t, time_discretization.gridInfo(time_stage_before_lift[i]).dt);
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(time_stage_before_lift[i]).dt+time_discretization.gridInfoLift(i).dt, dt);
//     EXPECT_FALSE(time_discretization.isSTOEnabledLift(i));
//   }
//   std::vector<int> time_stage_before_events;
//   for (auto e :  time_stage_before_impulse) {
//     time_stage_before_events.push_back(e);
//   }
//   for (auto e :  time_stage_before_lift) {
//     time_stage_before_events.push_back(e);
//   }
//   std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
//   time_stage_before_events.push_back(N+1);
//   int contact_phase_ref = 0;
//   for (int i=0; i<=N; ++i) {
//     EXPECT_EQ(time_discretization.contactPhase(i), contact_phase_ref);
//     EXPECT_EQ(time_discretization.gridInfo(i).phase, contact_phase_ref);
//     if (i == time_stage_before_events[contact_phase_ref]) {
//       ++contact_phase_ref;
//     }
//   }
//   for (int i=0; i<time_discretization.N_impulse(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeImpulse(i)+1, 
//               time_discretization.timeStageAfterImpulse(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterImpulse(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterImpulse(i)));
//   }
//   for (int i=0; i<time_discretization.N_lift(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeLift(i)+1, 
//               time_discretization.timeStageAfterLift(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterLift(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterLift(i)));
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeImpulse(time_discretization.impulseIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.impulseIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeLift(time_discretization.liftIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.liftIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<=N; ++i) {
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(i).t, t+i*dt);
//   }
//   for (int i=0; i<N; ++i) {
//     if (!time_discretization.isTimeStageBeforeImpulse(i) && !time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_DOUBLE_EQ(time_discretization.gridInfo(i).dt, dt);
//     }
//   }
//   const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
//   int impulse_index = 0;
//   int lift_index = 0;
//   for (int event_index=0; event_index<num_events; ++event_index) {
//     EXPECT_FALSE(time_discretization.eventType(event_index)==DiscreteEventType::None);
//     if (time_discretization.eventType(event_index) == DiscreteEventType::Impact) {
//       EXPECT_EQ(time_discretization.eventIndexImpulse(impulse_index), event_index);
//       ++impulse_index;
//     }
//     else {
//       EXPECT_EQ(time_discretization.eventIndexLift(lift_index), event_index);
//       ++lift_index;
//     }
//   }
//   EXPECT_NO_THROW(
//     std::cout << time_discretization << std::endl;
//   );
// }


// TEST_P(TimeDiscretizationTest, discretizeGridBased_switchingTimesOnGrids) {
//   TimeDiscretization time_discretization(T, N, max_num_events);
//   const auto robot = GetParam();
//   const auto contact_sequence = createContactSequenceOnGrid(robot);
//   time_discretization.discretize(contact_sequence, t);
//   EXPECT_EQ(time_discretization.N(), N-max_num_events);
//   EXPECT_EQ(time_discretization.N_impulse(), contact_sequence->numImpulseEvents());
//   EXPECT_EQ(time_discretization.N_lift(), contact_sequence->numLiftEvents());
//   double ti = t;
//   for (int i=0; i<time_discretization.N(); ++i) {
//     EXPECT_NEAR(time_discretization.gridInfo(i).dt, dt, min_dt);
//     EXPECT_NEAR(time_discretization.gridInfo(i).t, ti, min_dt);
//     ti += dt;
//     if (time_discretization.isTimeStageBeforeImpulse(i) || time_discretization.isTimeStageBeforeLift(i)) {
//       ti += dt;
//     }
//   }
//   EXPECT_DOUBLE_EQ(time_discretization.gridInfo(time_discretization.N()).t, t+T);
//   EXPECT_NO_THROW(
//     std::cout << time_discretization << std::endl;
//   );
// }


// TEST_P(TimeDiscretizationTest, discretizeGridBased_eventTimesAreLargerThanHorizon) {
//   const double T_short = 0.5 * T;
//   const double dt_short = T_short / N;
//   TimeDiscretization time_discretization(T_short, N, max_num_events);
//   const auto robot = GetParam();
//   const auto contact_sequence = createContactSequence(robot);
//   time_discretization.discretize(contact_sequence, t);
//   EXPECT_EQ(time_discretization.N(), N);
//   EXPECT_TRUE(time_discretization.N_impulse() <= contact_sequence->numImpulseEvents());
//   EXPECT_TRUE(time_discretization.N_lift() <= contact_sequence->numLiftEvents());
//   std::vector<double> t_impulse, t_lift;
//   std::vector<int> time_stage_before_impulse, time_stage_before_lift;
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     if (contact_sequence->impulseTime(i) > t+T_short-min_dt) {
//       break;
//     }
//     t_impulse.push_back(contact_sequence->impulseTime(i));
//     time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt_short));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     if (contact_sequence->liftTime(i) > t+T_short-min_dt) {
//       break;
//     }
//     t_lift.push_back(contact_sequence->liftTime(i));
//     time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt_short));
//   }
//   const int N_impulse = t_impulse.size();
//   const int N_lift = t_lift.size();
//   EXPECT_EQ(time_discretization.N_impulse(), N_impulse);
//   EXPECT_EQ(time_discretization.N_lift(), N_lift);
//   for (int i=0; i<N_impulse; ++i) {
//     EXPECT_EQ(time_stage_before_impulse[i], time_discretization.timeStageBeforeImpulse(i));
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.gridInfoImpulse(i).t);
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.gridInfoAux(i).t);
//     EXPECT_DOUBLE_EQ(t_impulse[i]-time_stage_before_impulse[i]*dt_short-t, 
//                      time_discretization.gridInfo(time_stage_before_impulse[i]).dt);
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(time_stage_before_impulse[i]).dt
//                      +time_discretization.gridInfoAux(i).dt, dt_short);
//     EXPECT_FALSE(time_discretization.isSTOEnabledImpulse(i));
//   }
//   for (int i=0; i<N_lift; ++i) {
//     EXPECT_EQ(time_stage_before_lift[i], time_discretization.timeStageBeforeLift(i));
//     EXPECT_DOUBLE_EQ(t_lift[i], time_discretization.gridInfoLift(i).t);
//     EXPECT_DOUBLE_EQ(t_lift[i]-time_stage_before_lift[i]*dt_short-t, 
//                      time_discretization.gridInfo(time_stage_before_lift[i]).dt);
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(time_stage_before_lift[i]).dt
//                      +time_discretization.gridInfoLift(i).dt, dt_short);
//     EXPECT_FALSE(time_discretization.isSTOEnabledLift(i));
//   }
//   std::vector<int> time_stage_before_events;
//   for (auto e :  time_stage_before_impulse) {
//     time_stage_before_events.push_back(e);
//   }
//   for (auto e :  time_stage_before_lift) {
//     time_stage_before_events.push_back(e);
//   }
//   std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
//   time_stage_before_events.push_back(N+1);
//   int contact_phase_ref = 0;
//   for (int i=0; i<=N; ++i) {
//     EXPECT_EQ(time_discretization.contactPhase(i), contact_phase_ref);
//     EXPECT_EQ(time_discretization.gridInfo(i).phase, contact_phase_ref);
//     if (i == time_stage_before_events[contact_phase_ref]) {
//       ++contact_phase_ref;
//     }
//   }
//   for (int i=0; i<time_discretization.N_impulse(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeImpulse(i)+1, 
//               time_discretization.timeStageAfterImpulse(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterImpulse(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterImpulse(i)));
//   }
//   for (int i=0; i<time_discretization.N_lift(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeLift(i)+1, 
//               time_discretization.timeStageAfterLift(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterLift(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterLift(i)));
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeImpulse(time_discretization.impulseIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.impulseIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeLift(time_discretization.liftIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.liftIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<=N; ++i) {
//     EXPECT_DOUBLE_EQ(time_discretization.gridInfo(i).t, t+i*dt_short);
//   }
//   for (int i=0; i<N; ++i) {
//     if (!time_discretization.isTimeStageBeforeImpulse(i) && !time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_DOUBLE_EQ(time_discretization.gridInfo(i).dt, dt_short);
//     }
//   }
//   const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
//   int impulse_index = 0;
//   int lift_index = 0;
//   for (int event_index=0; event_index<num_events; ++event_index) {
//     EXPECT_FALSE(time_discretization.eventType(event_index)==DiscreteEventType::None);
//     if (time_discretization.eventType(event_index) == DiscreteEventType::Impact) {
//       EXPECT_EQ(time_discretization.eventIndexImpulse(impulse_index), event_index);
//       ++impulse_index;
//     }
//     else {
//       EXPECT_EQ(time_discretization.eventIndexLift(lift_index), event_index);
//       ++lift_index;
//     }
//   }
//   EXPECT_NO_THROW(
//     std::cout << time_discretization << std::endl;
//   );
//   const auto time_steps = time_discretization.timeSteps();
//   int j=0;
//   for (int i=0; i<time_discretization.N(); ++i, ++j) {
//     EXPECT_DOUBLE_EQ(time_steps[j], time_discretization.gridInfo(i).dt);
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       ++j;
//       EXPECT_DOUBLE_EQ(time_steps[j], 
//                        time_discretization.gridInfoAux(time_discretization.impulseIndexAfterTimeStage(i)).dt);
//     }
//     else if (time_discretization.isTimeStageBeforeLift(i)) {
//       ++j;
//       EXPECT_DOUBLE_EQ(time_steps[j], 
//                        time_discretization.gridInfoLift(time_discretization.liftIndexAfterTimeStage(i)).dt);
//     }
//   }
// }


// TEST_P(TimeDiscretizationTest, discretizePhaseBased) {
//   TimeDiscretization time_discretization(T, N, max_num_events);
//   time_discretization.setDiscretizationMethod(DiscretizationMethod::PhaseBased);
//   const auto robot = GetParam();
//   const auto contact_sequence = createContactSequence(robot);
//   time_discretization.discretize(contact_sequence, t);
//   EXPECT_EQ(time_discretization.N(), N);
//   EXPECT_EQ(time_discretization.N_impulse(), contact_sequence->numImpulseEvents());
//   EXPECT_EQ(time_discretization.N_lift(), contact_sequence->numLiftEvents());
//   std::vector<double> t_impulse, t_lift;
//   std::vector<int> time_stage_before_impulse, time_stage_before_lift;
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     t_impulse.push_back(contact_sequence->impulseTime(i));
//     time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     t_lift.push_back(contact_sequence->liftTime(i));
//     time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt));
//   }
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     EXPECT_EQ(time_stage_before_impulse[i], time_discretization.timeStageBeforeImpulse(i));
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.gridInfoImpulse(i).t);
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.impulseTime(i));
//     EXPECT_FALSE(time_discretization.isSTOEnabledImpulse(i));
//     EXPECT_EQ(time_discretization.isSTOEnabledImpulse(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexImpulse(i)));
//     EXPECT_EQ(time_discretization.isSTOEnabledImpulse(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexImpulse(i)));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     EXPECT_EQ(time_stage_before_lift[i], time_discretization.timeStageBeforeLift(i));
//     EXPECT_DOUBLE_EQ(t_lift[i], time_discretization.gridInfoLift(i).t);
//     EXPECT_DOUBLE_EQ(t_lift[i], time_discretization.liftTime(i));
//     EXPECT_FALSE(time_discretization.isSTOEnabledLift(i));
//     EXPECT_EQ(time_discretization.isSTOEnabledLift(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexLift(i)));
//   }
//   EXPECT_EQ(time_discretization.isSTOEnabledPhase(0), 
//             time_discretization.isSTOEnabledEvent(0));
//   for (int i=1; i<contact_sequence->numDiscreteEvents(); ++i) {
//     EXPECT_EQ(time_discretization.isSTOEnabledPhase(i), 
//               (time_discretization.isSTOEnabledEvent(i-1)||time_discretization.isSTOEnabledEvent(i)));
//   }
//   EXPECT_EQ(time_discretization.isSTOEnabledPhase(contact_sequence->numDiscreteEvents()), 
//             time_discretization.isSTOEnabledEvent(contact_sequence->numDiscreteEvents()-1));
//   std::vector<int> time_stage_before_events;
//   for (auto e :  time_stage_before_impulse) {
//     time_stage_before_events.push_back(e);
//   }
//   for (auto e :  time_stage_before_lift) {
//     time_stage_before_events.push_back(e);
//   }
//   std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
//   time_stage_before_events.push_back(N+1);
//   int contact_phase_ref = 0;
//   for (int i=0; i<=N; ++i) {
//     EXPECT_EQ(time_discretization.contactPhase(i), contact_phase_ref);
//     EXPECT_EQ(time_discretization.gridInfo(i).phase, contact_phase_ref);
//     if (i == time_stage_before_events[contact_phase_ref]) {
//       ++contact_phase_ref;
//     }
//   }
//   for (int i=0; i<time_discretization.N_impulse(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeImpulse(i)+1, 
//               time_discretization.timeStageAfterImpulse(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterImpulse(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterImpulse(i)));
//   }
//   for (int i=0; i<time_discretization.N_lift(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeLift(i)+1, 
//               time_discretization.timeStageAfterLift(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterLift(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterLift(i)));
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeImpulse(time_discretization.impulseIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.impulseIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeLift(time_discretization.liftIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.liftIndexAfterTimeStage(i), -1);
//     }
//   }
//   const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
//   int impulse_index = 0;
//   int lift_index = 0;
//   int time_stage_before_event = 0;
//   double t_prev_event = t;
//   double dt_prev_aux = time_discretization.gridInfo(0).dt;
//   for (int event_index=0; event_index<num_events; ++event_index) {
//     ASSERT_FALSE(time_discretization.eventType(event_index)==DiscreteEventType::None);
//     if (time_discretization.eventType(event_index) == DiscreteEventType::Impact) {
//       EXPECT_EQ(time_discretization.eventIndexImpulse(impulse_index), event_index);
//       const int grids_phase = time_discretization.timeStageBeforeImpulse(impulse_index) 
//                               - time_stage_before_event + 1;
//       EXPECT_EQ(time_discretization.num_grids_in_phase(event_index), grids_phase);
//       const double dt_phase = (time_discretization.impulseTime(impulse_index)-t_prev_event) / grids_phase;
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeImpulse(impulse_index); ++stage) {
//         EXPECT_DOUBLE_EQ(time_discretization.gridInfo(stage).dt, dt_phase);
//       }
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeImpulse(impulse_index)-1; ++stage) {
//         EXPECT_NEAR(time_discretization.gridInfo(stage+1).t-time_discretization.gridInfo(stage).t, 
//                     dt_phase, min_dt);
//       }
//       EXPECT_NEAR(time_discretization.impulseTime(impulse_index)-time_discretization.gridInfo(time_discretization.timeStageBeforeImpulse(impulse_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_NEAR(time_discretization.gridInfoImpulse(impulse_index).t-time_discretization.gridInfo(time_discretization.timeStageBeforeImpulse(impulse_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
//       time_stage_before_event = time_discretization.timeStageBeforeImpulse(impulse_index);
//       t_prev_event = time_discretization.gridInfoImpulse(impulse_index).t;
//       dt_prev_aux = time_discretization.gridInfoAux(impulse_index).dt;
//       ++impulse_index;
//     }
//     else {
//       EXPECT_EQ(time_discretization.eventIndexLift(lift_index), event_index);
//       const int grids_phase = time_discretization.timeStageBeforeLift(lift_index) 
//                               - time_stage_before_event + 1;
//       EXPECT_EQ(time_discretization.num_grids_in_phase(event_index), grids_phase);
//       const double dt_phase = (time_discretization.liftTime(lift_index)-t_prev_event) / grids_phase;
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeLift(lift_index); ++stage) {
//         EXPECT_DOUBLE_EQ(time_discretization.gridInfo(stage).dt, dt_phase);
//       }
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeLift(lift_index)-1; ++stage) {
//         EXPECT_NEAR(time_discretization.gridInfo(stage+1).t-time_discretization.gridInfo(stage).t, dt_phase, min_dt);
//       }
//       EXPECT_NEAR(time_discretization.liftTime(lift_index)-time_discretization.gridInfo(time_discretization.timeStageBeforeLift(lift_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_NEAR(time_discretization.gridInfoLift(lift_index).t-time_discretization.gridInfo(time_discretization.timeStageBeforeLift(lift_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
//       time_stage_before_event = time_discretization.timeStageBeforeLift(lift_index);
//       t_prev_event = time_discretization.liftTime(lift_index);
//       dt_prev_aux = time_discretization.gridInfoLift(lift_index).dt;
//       ++lift_index;
//     }
//   }
//   EXPECT_NO_THROW(
//     std::cout << time_discretization << std::endl;
//   );
// }


// TEST_P(TimeDiscretizationTest, discretizePhaseBased_eventTimesAreLargerThanHorizon) {
//   const double T_short = 0.5 * T;
//   const double dt_short = T_short / N;
//   TimeDiscretization time_discretization(T_short, N, max_num_events);
//   time_discretization.setDiscretizationMethod(DiscretizationMethod::PhaseBased);
//   const auto robot = GetParam();
//   const auto contact_sequence = createContactSequence(robot);
//   time_discretization.discretize(contact_sequence, t);
//   EXPECT_EQ(time_discretization.N(), N);
//   std::vector<double> t_impulse, t_lift;
//   std::vector<int> time_stage_before_impulse, time_stage_before_lift;
//   for (int i=0; i<contact_sequence->numImpulseEvents(); ++i) {
//     if (contact_sequence->impulseTime(i) >= t+T_short-min_dt) {
//       break;
//     }
//     t_impulse.push_back(contact_sequence->impulseTime(i));
//     time_stage_before_impulse.push_back(std::floor((t_impulse[i]-t)/dt_short));
//   }
//   for (int i=0; i<contact_sequence->numLiftEvents(); ++i) {
//     if (contact_sequence->liftTime(i) >= t+T_short-min_dt) {
//       break;
//     }
//     t_lift.push_back(contact_sequence->liftTime(i));
//     time_stage_before_lift.push_back(std::floor((t_lift[i]-t)/dt_short));
//   }
//   const int N_impulse = t_impulse.size();
//   const int N_lift = t_lift.size();
//   EXPECT_EQ(time_discretization.N_impulse(), N_impulse);
//   EXPECT_EQ(time_discretization.N_lift(), N_lift);
//   for (int i=0; i<N_impulse; ++i) {
//     EXPECT_EQ(time_stage_before_impulse[i], time_discretization.timeStageBeforeImpulse(i));
//     EXPECT_DOUBLE_EQ(t_impulse[i], time_discretization.impulseTime(i));
//     EXPECT_FALSE(time_discretization.isSTOEnabledImpulse(i));
//     EXPECT_EQ(time_discretization.isSTOEnabledImpulse(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexImpulse(i)));
//     EXPECT_EQ(time_discretization.isSTOEnabledImpulse(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexImpulse(i)));
//   }
//   for (int i=0; i<N_lift; ++i) {
//     EXPECT_EQ(time_stage_before_lift[i], time_discretization.timeStageBeforeLift(i));
//     EXPECT_DOUBLE_EQ(t_lift[i], time_discretization.liftTime(i));
//     EXPECT_FALSE(time_discretization.isSTOEnabledLift(i));
//     EXPECT_EQ(time_discretization.isSTOEnabledLift(i), 
//               time_discretization.isSTOEnabledEvent(time_discretization.eventIndexLift(i)));
//   }
//   EXPECT_EQ(time_discretization.isSTOEnabledPhase(0), 
//             time_discretization.isSTOEnabledEvent(0));
//   for (int i=1; i<N_lift+N_impulse; ++i) {
//     EXPECT_EQ(time_discretization.isSTOEnabledPhase(i), 
//               (time_discretization.isSTOEnabledEvent(i-1)||time_discretization.isSTOEnabledEvent(i)));
//   }
//   EXPECT_EQ(time_discretization.isSTOEnabledPhase(N_lift+N_impulse), 
//             time_discretization.isSTOEnabledEvent(N_lift+N_impulse-1));
//   std::vector<int> time_stage_before_events;
//   for (auto e :  time_stage_before_impulse) {
//     time_stage_before_events.push_back(e);
//   }
//   for (auto e :  time_stage_before_lift) {
//     time_stage_before_events.push_back(e);
//   }
//   std::sort(time_stage_before_events.begin(), time_stage_before_events.end());
//   time_stage_before_events.push_back(N+1);
//   int contact_phase_ref = 0;
//   for (int i=0; i<=N; ++i) {
//     EXPECT_EQ(time_discretization.contactPhase(i), contact_phase_ref);
//     if (i == time_stage_before_events[contact_phase_ref]) {
//       ++contact_phase_ref;
//     }
//   }
//   for (int i=0; i<time_discretization.N_impulse(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeImpulse(i)+1, 
//               time_discretization.timeStageAfterImpulse(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterImpulse(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterImpulse(i)));
//   }
//   for (int i=0; i<time_discretization.N_lift(); ++i) {
//     EXPECT_EQ(time_discretization.timeStageBeforeLift(i)+1, 
//               time_discretization.timeStageAfterLift(i));
//     EXPECT_EQ(time_discretization.contactPhaseAfterLift(i), 
//               time_discretization.contactPhase(time_discretization.timeStageAfterLift(i)));
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeImpulse(time_discretization.impulseIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.impulseIndexAfterTimeStage(i), -1);
//     }
//   }
//   for (int i=0; i<N; ++i) {
//     if (time_discretization.isTimeStageBeforeLift(i)) {
//       EXPECT_EQ(time_discretization.timeStageBeforeLift(time_discretization.liftIndexAfterTimeStage(i)), i);
//     }
//     else {
//       EXPECT_EQ(time_discretization.liftIndexAfterTimeStage(i), -1);
//     }
//   }
//   const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
//   int impulse_index = 0;
//   int lift_index = 0;
//   int time_stage_before_event = 0;
//   double t_prev_event = t;
//   double dt_prev_aux = time_discretization.gridInfo(0).dt;
//   for (int event_index=0; event_index<num_events; ++event_index) {
//     ASSERT_FALSE(time_discretization.eventType(event_index)==DiscreteEventType::None);
//     if (time_discretization.eventType(event_index) == DiscreteEventType::Impact) {
//       EXPECT_EQ(time_discretization.eventIndexImpulse(impulse_index), event_index);
//       const int grids_phase = time_discretization.timeStageBeforeImpulse(impulse_index) 
//                               - time_stage_before_event + 1;
//       EXPECT_EQ(time_discretization.num_grids_in_phase(event_index), grids_phase);
//       const double dt_phase = (time_discretization.impulseTime(impulse_index)-t_prev_event) / grids_phase;
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeImpulse(impulse_index); ++stage) {
//         EXPECT_DOUBLE_EQ(time_discretization.gridInfo(stage).dt, dt_phase);
//       }
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeImpulse(impulse_index)-1; ++stage) {
//         EXPECT_NEAR(time_discretization.gridInfo(stage+1).t-time_discretization.gridInfo(stage).t, dt_phase, min_dt);
//       }
//       EXPECT_NEAR(time_discretization.impulseTime(impulse_index)-time_discretization.gridInfo(time_discretization.timeStageBeforeImpulse(impulse_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
//       time_stage_before_event = time_discretization.timeStageBeforeImpulse(impulse_index);
//       t_prev_event = time_discretization.impulseTime(impulse_index);
//       dt_prev_aux = time_discretization.gridInfoAux(impulse_index).dt;
//       ++impulse_index;
//     }
//     else {
//       EXPECT_EQ(time_discretization.eventIndexLift(lift_index), event_index);
//       const int grids_phase = time_discretization.timeStageBeforeLift(lift_index) 
//                               - time_stage_before_event + 1;
//       EXPECT_EQ(time_discretization.num_grids_in_phase(event_index), grids_phase);
//       const double dt_phase = (time_discretization.liftTime(lift_index)-t_prev_event) / grids_phase;
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeLift(lift_index); ++stage) {
//         EXPECT_DOUBLE_EQ(time_discretization.gridInfo(stage).dt, dt_phase);
//       }
//       for (int stage=time_stage_before_event+1; 
//             stage<time_discretization.timeStageBeforeLift(lift_index)-1; ++stage) {
//         EXPECT_NEAR(time_discretization.gridInfo(stage+1).t-time_discretization.gridInfo(stage).t, dt_phase, min_dt);
//       }
//       EXPECT_NEAR(time_discretization.liftTime(lift_index)-time_discretization.gridInfo(time_discretization.timeStageBeforeLift(lift_index)).t, 
//                   dt_phase, min_dt);
//       EXPECT_DOUBLE_EQ(dt_prev_aux, dt_phase);
//       time_stage_before_event = time_discretization.timeStageBeforeLift(lift_index);
//       t_prev_event = time_discretization.liftTime(lift_index);
//       dt_prev_aux = time_discretization.gridInfoLift(lift_index).dt;
//       ++lift_index;
//     }
//   }
//   EXPECT_NO_THROW(
//     std::cout << time_discretization << std::endl;
//   );
//   const auto time_steps = time_discretization.timeSteps();
//   int j=0;
//   for (int i=0; i<time_discretization.N(); ++i, ++j) {
//     EXPECT_DOUBLE_EQ(time_steps[j], time_discretization.gridInfo(i).dt);
//     if (time_discretization.isTimeStageBeforeImpulse(i)) {
//       ++j;
//       EXPECT_DOUBLE_EQ(time_steps[j], 
//                        time_discretization.gridInfoAux(time_discretization.impulseIndexAfterTimeStage(i)).dt);
//     }
//     else if (time_discretization.isTimeStageBeforeLift(i)) {
//       ++j;
//       EXPECT_DOUBLE_EQ(time_steps[j], 
//                        time_discretization.gridInfoLift(time_discretization.liftIndexAfterTimeStage(i)).dt);
//     }
//   }
// }


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, TimeDiscretizationTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateQuadrupedalRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}