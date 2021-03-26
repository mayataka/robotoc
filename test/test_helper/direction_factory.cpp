#include "direction_factory.hpp"

#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"


namespace idocp {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N, 
                          const int max_num_impulse) {
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    d[i].setRandom();
  }
  return d;
}


Direction CreateDirection(const Robot& robot, const ContactSequence& contact_sequence, 
                          const double T, const int N, const int max_num_impulse, 
                          const double t, const bool is_parnmpc) {
  if (robot.maxPointContacts() == 0) {
    return CreateDirection(robot, N, max_num_impulse);
  }
  else {
    if (is_parnmpc) {
      ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_impulse);
      parnmpc_discretizer.discretizeOCP(contact_sequence, t);
      Direction d(robot, N, max_num_impulse);
      for (int i=0; i<N; ++i) {
        d[i].setRandom(contact_sequence.contactStatus(parnmpc_discretizer.contactPhase(i)));
      }
      const int num_impulse = contact_sequence.numImpulseEvents();
      for (int i=0; i<num_impulse; ++i) {
        d.impulse[i].setRandom(contact_sequence.impulseStatus(i));
      }
      for (int i=0; i<num_impulse; ++i) {
        d.aux[i].setRandom(contact_sequence.contactStatus(parnmpc_discretizer.contactPhaseBeforeImpulse(i)));
        if (parnmpc_discretizer.timeStageBeforeImpulse(i) >= 0) {
          d.aux[i].setImpulseStatus(contact_sequence.impulseStatus(i));
          d.aux[i].dxi().setRandom();
        }
      }
      const int num_lift = contact_sequence.numLiftEvents();
      for (int i=0; i<num_lift; ++i) {
        d.lift[i].setRandom(contact_sequence.contactStatus(parnmpc_discretizer.contactPhaseBeforeLift(i)));
      }
      return d;
    }
    else {
      OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
      ocp_discretizer.discretizeOCP(contact_sequence, t);
      Direction d(robot, N, max_num_impulse);
      for (int i=0; i<=N; ++i) {
        d[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhase(i)));
      }
      const int num_impulse = contact_sequence.numImpulseEvents();
      for (int i=0; i<num_impulse; ++i) {
        d.impulse[i].setRandom(contact_sequence.impulseStatus(i));
      }
      for (int i=0; i<num_impulse; ++i) {
        d.aux[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterImpulse(i)));
      }
      for (int i=0; i<num_impulse; ++i) {
        const int time_stage_before_impulse = ocp_discretizer.timeStageBeforeImpulse(i);
        if (time_stage_before_impulse-1 >= 0) {
          d[time_stage_before_impulse-1].setImpulseStatus(contact_sequence.impulseStatus(i));
          d[time_stage_before_impulse-1].dxi().setRandom();
        }
      }
      const int num_lift = contact_sequence.numLiftEvents();
      for (int i=0; i<num_lift; ++i) {
        d.lift[i].setRandom(contact_sequence.contactStatus(ocp_discretizer.contactPhaseAfterLift(i)));
      }
      return d;
    }
  }
}

} // namespace testhelper
} // namespace idocp