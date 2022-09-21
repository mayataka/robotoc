#include "solution_factory.hpp"

#include "robotoc/hybrid/time_discretization.hpp"


namespace robotoc {
namespace testhelper {

Solution CreateSolution(const Robot& robot, const int N, 
                        const int max_num_impulse) {
  Solution s(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution CreateSolution(const Robot& robot,  
                        const std::shared_ptr<ContactSequence>& contact_sequence, 
                        const double T, const int N, 
                        const int max_num_impulse, const double t) {
  if (robot.maxNumContacts() == 0) {
    return CreateSolution(robot, N, max_num_impulse);
  }
  else {
    TimeDiscretization time_discretization(T, N, max_num_impulse);
    time_discretization.discretize(contact_sequence, t);
    Solution s(robot, N, max_num_impulse);
    for (int i=0; i<=N; ++i) {
      s[i].setRandom(robot, contact_sequence->contactStatus(time_discretization.contactPhase(i)));
    }
    const int num_impulse = contact_sequence->numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      s.impulse[i].setRandom(robot, contact_sequence->impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      s.aux[i].setRandom(robot, contact_sequence->contactStatus(time_discretization.contactPhaseAfterImpulse(i)));
    }
    for (int i=0; i<num_impulse; ++i) {
      const int time_stage_before_impulse = time_discretization.timeStageBeforeImpulse(i);
      if (time_stage_before_impulse-1 >= 0) {
        s[time_stage_before_impulse-1].setSwitchingConstraint(contact_sequence->impulseStatus(i));
        s[time_stage_before_impulse-1].xi_stack().setRandom();
      }
    }
    const int num_lift = contact_sequence->numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      s.lift[i].setRandom(robot, contact_sequence->contactStatus(time_discretization.contactPhaseAfterLift(i)));
    }
    return s;
  }
}

} // namespace testhelper
} // namespace robotoc