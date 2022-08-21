#include "direction_factory.hpp"

#include "robotoc/hybrid/time_discretization.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N, 
                          const int max_num_impulse) {
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    d[i].setRandom();
  }
  return d;
}


Direction CreateDirection(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence, 
                          const double T, const int N, 
                          const int max_num_impulse, const double t) {
  if (robot.maxNumContacts() == 0) {
    return CreateDirection(robot, N, max_num_impulse);
  }
  else {
    TimeDiscretization time_discretization(T, N, max_num_impulse);
    time_discretization.discretize(contact_sequence, t);
    Direction d(robot, N, max_num_impulse);
    for (int i=0; i<=N; ++i) {
      d[i].setRandom(contact_sequence->contactStatus(time_discretization.contactPhase(i)));
    }
    const int num_impulse = contact_sequence->numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      d.impulse[i].setRandom(contact_sequence->impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      d.aux[i].setRandom(contact_sequence->contactStatus(time_discretization.contactPhaseAfterImpulse(i)));
    }
    for (int i=0; i<num_impulse; ++i) {
      const int time_stage_before_impulse = time_discretization.timeStageBeforeImpulse(i);
      if (time_stage_before_impulse-1 >= 0) {
        d[time_stage_before_impulse-1].setImpulseStatus(contact_sequence->impulseStatus(i));
        d[time_stage_before_impulse-1].dxi().setRandom();
      }
    }
    const int num_lift = contact_sequence->numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      d.lift[i].setRandom(contact_sequence->contactStatus(time_discretization.contactPhaseAfterLift(i)));
    }
    return d;
  }
}

} // namespace testhelper
} // namespace robotoc