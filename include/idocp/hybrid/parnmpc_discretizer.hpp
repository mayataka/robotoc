#ifndef IDOCP_PARNMPC_DISCRETIZER_HPP_
#define IDOCP_PARNMPC_DISCRETIZER_HPP_

#include "idocp/hybrid/contact_sequence.hpp"

#include <vector>
#include <limits>

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

  int N_impulse() const;

  int N_lift() const; 

  int N_all() const;

  int N_ideal() const;

  int contactPhase(const int time_stage) const;

  int impulseIndexBeforeTimeStage(const int time_stage) const;

  int impulseIndexAfterTimeStage(const int time_stage) const;

  int liftIndexBeforeTimeStage(const int time_stage) const;

  int liftIndexAfterTimeStage(const int time_stage) const;

  int timeStageBeforeImpulse(const int impulse_index) const;

  int timeStageAfterImpulse(const int impulse_index) const;

  int timeStageBeforeLift(const int lift_index) const;

  int timeStageAfterLift(const int lift_index) const;

  int contactPhaseBeforeImpulse(const int impulse_index) const;

  int contactPhaseBeforeLift(const int lift_index) const;

  bool isTimeStageBeforeImpulse(const int time_stage) const;

  bool isTimeStageAfterImpulse(const int time_stage) const;

  bool isTimeStageBeforeLift(const int time_stage) const;

  bool isTimeStageAfterLift(const int time_stage) const;

  double t(const int time_stage) const;

  double t_impulse(const int impulse_index) const;

  double t_lift(const int lift_index) const;

  double dt(const int time_stage) const;

  double dt_aux(const int impulse_index) const;

  double dt_lift(const int lift_index) const;

  bool isWellDefined() const;

private:
  double T_, dt_ideal_, max_dt_;
  int N_, N_ideal_, N_impulse_, N_lift_, max_events_;
  std::vector<int> contact_phase_index_from_time_stage_, 
                   impulse_index_before_time_stage_, 
                   lift_index_before_time_stage_, time_stage_after_impulse_, 
                   time_stage_after_lift_;
  std::vector<bool> is_time_stage_after_impulse_, is_time_stage_after_lift_;
  std::vector<double> t_, t_impulse_, t_lift_, dt_, dt_aux_, dt_lift_;

  static constexpr double min_dt_ 
      = std::sqrt(std::numeric_limits<double>::epsilon());

  void countDiscreteEvents(const ContactSequence& contact_sequence, 
                           const double t);

  void countTimeSteps(const double t);

  void countTimeStages();

  void countContactPhase();

};

} // namespace idocp

#include "idocp/hybrid/parnmpc_discretizer.hxx"

#endif // IDOCP_PARNMPC_DISCRETIZER_HPP_ 