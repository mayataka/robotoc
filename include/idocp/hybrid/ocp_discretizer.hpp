#ifndef IDOCP_OCP_DISCRETIZER_HPP_
#define IDOCP_OCP_DISCRETIZER_HPP_

#include "idocp/hybrid/contact_sequence.hpp"

#include <vector>
#include <limits>

namespace idocp {

class OCPDiscretizer {
public:
  OCPDiscretizer(const double T, const int N, const int max_events);

  ///
  /// @brief Default constructor. 
  ///
  OCPDiscretizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPDiscretizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPDiscretizer(const OCPDiscretizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPDiscretizer& operator=(const OCPDiscretizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPDiscretizer(OCPDiscretizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPDiscretizer& operator=(OCPDiscretizer&&) noexcept = default;

  void discretizeOCP(const ContactSequence& contact_sequence, const double t);

  int N() const;

  int numImpulseStages() const;

  int numLiftStages() const; 

  int contactPhase(const int time_stage) const;

  int impulseIndex(const int time_stage_before_impulse) const;

  int liftIndex(const int time_stage_before_lift) const;

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

  bool existImpulse() const;

  double t(const int time_stage) const;

  double t_impulse(const int impulse_index) const;

  double t_lift(const int lift_index) const;

  double dtau(const int time_stage) const;

  double dtau_aux(const int impulse_index) const;

  double dtau_lift(const int lift_index) const;

private:
  double T_;
  int N_, max_events_, num_impulse_stages_, num_lift_stages_;
  std::vector<int> contact_phase_index_from_time_stage_, 
                   impulse_index_from_time_stage_, lift_index_from_time_stage_, 
                   time_stage_before_impulse_, time_stage_before_lift_;
  std::vector<bool> is_time_stage_before_impulse_, is_time_stage_before_lift_;
  std::vector<double> t_, t_impulse_, t_lift_, dtau_, dtau_aux_, dtau_lift_;

  static constexpr double kMinDiscretizationSize
      = std::sqrt(std::numeric_limits<double>::epsilon());

  void countImpulseEvents(const ContactSequence& contact_sequence, 
                          const double t);

  void countLiftEvents(const ContactSequence& contact_sequence, const double t);

  void countContactPhase(const ContactSequence& contact_sequence);

  void countIsTimeStageBeforeEvents(const ContactSequence& contact_sequence);

  void countTime(const ContactSequence& contact_sequence, const double t);

};

} // namespace idocp

#include "idocp/hybrid/ocp_discretizer.hxx"

#endif // IDOCP_OCP_DISCRETIZER_HPP_ 