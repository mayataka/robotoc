#include "solution_factory.hpp"

#include "idocp/hybrid/ocp_discretizer.hpp"


namespace idocp {
namespace testhelper {

Solution CreateSolution(const Robot& robot, const int N, 
                        const int max_num_impulse) {
  Solution s(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    s[i].setRandom(robot);
  }
  return s;
}


Solution CreateSolution(const Robot& robot, const ContactSequence& contact_sequence, 
                        const double T, const int N, const int max_num_impulse, const double t) {
  if (robot.maxPointContacts() == 0) {
    return CreateSolution(robot, N, max_num_impulse);
  }
  else {
    OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
    ocp_discretizer.discretizeOCP(contact_sequence, t);
    Solution s(robot, N, max_num_impulse);
    for (int i=0; i<=N; ++i) {
      s[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)));
    }
    const int num_impulse = contact_sequence.numImpulseEvents();
    for (int i=0; i<num_impulse; ++i) {
      s.impulse[i].setRandom(robot, contact_sequence.impulseStatus(i));
    }
    for (int i=0; i<num_impulse; ++i) {
      s.aux[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterImpulse(i)));
    }
    for (int i=0; i<num_impulse; ++i) {
      const int time_stage_before_impulse = ocp_discretizer.timeStageBeforeImpulse(i);
      if (time_stage_before_impulse-1 >= 0) {
        s[time_stage_before_impulse-1].setImpulseStatus(contact_sequence.impulseStatus(i));
        s[time_stage_before_impulse-1].xi_stack().setRandom();
      }
    }
    const int num_lift = contact_sequence.numLiftEvents();
    for (int i=0; i<num_lift; ++i) {
      s.lift[i].setRandom(robot, contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterLift(i)));
    }
    return s;
  }
}

} // namespace testhelper
} // namespace idocp