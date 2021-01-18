#ifndef IDOCP_PARNMPC_DISCRETIZER_HPP_
#define IDOCP_PARNMPC_DISCRETIZER_HPP_

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"


namespace idocp {

class ParNMPCDiscretizer {
public:
  ParNMPCDiscretizer(const double T, const int N, const int max_events);

  ///
  /// @brief Default constructor. 
  ///
  ParNMPCDiscretizer();

  ///
  /// @brief Destructor. 
  ///
  ~ParNMPCDiscretizer();

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCDiscretizer(const ParNMPCDiscretizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCDiscretizer& operator=(const ParNMPCDiscretizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCDiscretizer(ParNMPCDiscretizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCDiscretizer& operator=(ParNMPCDiscretizer&&) noexcept = default;

  void discretizeOCP(const ContactSequence& contact_sequence, const double t);

  int N() const;

  int numImpulseStages() const;

  int numLiftStages() const; 

  int contactPhase(const int time_stage) const;

  int impulseIndexAfterTimeStage(const int time_stage) const;

  int liftIndexAfterTimeStage(const int time_stage) const;

  int timeStageBeforeImpulse(const int impulse_index) const;

  int timeStageAfterImpulse(const int impulse_index) const;

  int timeStageBeforeLift(const int lift_index) const;

  int timeStageAfterLift(const int lift_index) const;

  int contactPhaseBeforeImpulse(const int impulse_index) const;

  int contactPhaseAfterImpulse(const int impulse_index) const;

  int contactPhaseBeforeLift(const int lift_index) const;

  int contactPhaseAfterLift(const int lift_index) const;

  bool isTimeStageBeforeImpulse(const int time_stage) const;

  bool isTimeStageAfterImpulse(const int time_stage) const;

  bool isTimeStageBeforeLift(const int time_stage) const;

  bool isTimeStageAfterLift(const int time_stage) const;

  double t(const int time_stage) const;

  double t_impulse(const int impulse_index) const;

  double t_lift(const int lift_index) const;

  double dtau(const int time_stage) const;

  double dtau_aux(const int impulse_index) const;

  double dtau_lift(const int lift_index) const;

private:
  OCPDiscretizer ocp_discretizer_;

};

} // namespace idocp

#include "idocp/hybrid/parnmpc_discretizer.hxx"

#endif // IDOCP_PARNMPC_DISCRETIZER_HPP_ 