#ifndef IDOCP_PARNMPC_DISCRETIZER_HXX_
#define IDOCP_PARNMPC_DISCRETIZER_HXX_

#include "idocp/hybrid/parnmpc_discretizer.hpp"


namespace idocp {

inline ParNMPCDiscretizer::ParNMPCDiscretizer(const double T, const int N, 
                                              const int max_events) 
  : ocp_discretizer_(T, N, max_events) {
}


inline ParNMPCDiscretizer::ParNMPCDiscretizer()
  : ocp_discretizer_() {
}


inline ParNMPCDiscretizer::~ParNMPCDiscretizer() {
}


inline void ParNMPCDiscretizer::discretizeOCP(
    const ContactSequence& contact_sequence, const double t) {
  ocp_discretizer_.discretizeOCP(contact_sequence, t);
}


inline int ParNMPCDiscretizer::N() const {
  return ocp_discretizer_.N();
}


inline int ParNMPCDiscretizer::numImpulseStages() const {
  return ocp_discretizer_.numImpulseStages();
}


inline int ParNMPCDiscretizer::numLiftStages() const {
  return ocp_discretizer_.numLiftStages();
}


inline int ParNMPCDiscretizer::contactPhase(const int time_stage) const {
  return ocp_discretizer_.contactPhase(time_stage+1);
}


inline int ParNMPCDiscretizer::impulseIndexAfterTimeStage(
    const int time_stage) const {
  return ocp_discretizer_.impulseIndexAfterTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::liftIndexAfterTimeStage(
    const int time_stage) const {
  return ocp_discretizer_.liftIndexAfterTimeStage(time_stage+1);
}


inline int ParNMPCDiscretizer::timeStageBeforeImpulse(
    const int impulse_index) const {
  return ocp_discretizer_.timeStageBeforeImpulse(impulse_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterImpulse(
    const int impulse_index) const {
  return ocp_discretizer_.timeStageAfterImpulse(impulse_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageBeforeLift(const int lift_index) const {
  return ocp_discretizer_.timeStageBeforeLift(lift_index) - 1;
}


inline int ParNMPCDiscretizer::timeStageAfterLift(const int lift_index) const {
  return ocp_discretizer_.timeStageAfterLift(lift_index) - 1;
}


inline int ParNMPCDiscretizer::contactPhaseBeforeImpulse(
    const int impulse_index) const {
  return ocp_discretizer_.contactPhaseBeforeImpulse(impulse_index);
}


inline int ParNMPCDiscretizer::contactPhaseAfterImpulse(
    const int impulse_index) const {
  return ocp_discretizer_.contactPhaseAfterImpulse(impulse_index);
}


inline int ParNMPCDiscretizer::contactPhaseBeforeLift(
    const int lift_index) const {
  return ocp_discretizer_.contactPhaseBeforeLift(lift_index);
}


inline int ParNMPCDiscretizer::contactPhaseAfterLift(
    const int lift_index) const {
  return ocp_discretizer_.contactPhaseAfterLift(lift_index);
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeImpulse(
    const int time_stage) const {
  return ocp_discretizer_.isTimeStageBeforeImpulse(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterImpulse(
    const int time_stage) const {
  return ocp_discretizer_.isTimeStageAfterImpulse(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageBeforeLift(
    const int time_stage) const {
  return ocp_discretizer_.isTimeStageBeforeLift(time_stage+1);
}


inline bool ParNMPCDiscretizer::isTimeStageAfterLift(
    const int time_stage) const {
  return ocp_discretizer_.isTimeStageAfterLift(time_stage+1);
}


inline double ParNMPCDiscretizer::t(const int time_stage) const {
  return ocp_discretizer_.t(time_stage+1);
}


inline double ParNMPCDiscretizer::t_impulse(const int impulse_index) const {
  return ocp_discretizer_.t_impulse(impulse_index);
}


inline double ParNMPCDiscretizer::t_lift(const int lift_index) const {
  return ocp_discretizer_.t_lift(lift_index);
}


inline double ParNMPCDiscretizer::dtau(const int time_stage) const {
  return ocp_discretizer_.dtau(time_stage+1);
}


inline double ParNMPCDiscretizer::dtau_aux(const int impulse_index) const {
  return ocp_discretizer_.dtau_aux(impulse_index);
}


inline double ParNMPCDiscretizer::dtau_lift(const int lift_index) const {
  return ocp_discretizer_.dtau_lift(lift_index);
}

} // namespace idocp

#endif // IDOCP_PARNMPC_DISCRETIZER_HXX_ 