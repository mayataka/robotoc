#ifndef ROBOTOC_HYBRID_OCP_DISCRETIZATION_HPP_ 
#define ROBOTOC_HYBRID_OCP_DISCRETIZATION_HPP_

#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <iostream>

#include "robotoc/hybrid/contact_sequence.hpp"


namespace robotoc {

///
/// @class HybridOCPDiscretization
/// @brief Non-uniform time discretization of the hybrid optimal control problem.
///
class HybridOCPDiscretization {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] T Length of the horizon.
  /// @param[in] N Ideal number of the discretization grids of the horizon. 
  /// @param[in] max_events Maximum number of each discrete events 
  /// (impulse and lift). 
  /// @note The actual number of the grids can differ from N depending on the 
  /// discrete events.
  ///
  HybridOCPDiscretization(const double T, const int N, const int max_events);

  ///
  /// @brief Default constructor. 
  ///
  HybridOCPDiscretization();

  ///
  /// @brief Destructor. 
  ///
  ~HybridOCPDiscretization();

  ///
  /// @brief Default copy constructor. 
  ///
  HybridOCPDiscretization(const HybridOCPDiscretization&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  HybridOCPDiscretization& operator=(const HybridOCPDiscretization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  HybridOCPDiscretization(HybridOCPDiscretization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  HybridOCPDiscretization& operator=(HybridOCPDiscretization&&) noexcept = default;

  ///
  /// @brief Discretizes the finite horizon taking into account the discrete 
  /// events. 
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  /// @param[in] t Initial time of the horizon.
  ///
  void discretize(const std::shared_ptr<ContactSequence>& contact_sequence, 
                  const double t);

  ///
  /// @return Number of the time stages on the horizon. 
  ///
  int N() const;

  ///
  /// @return Number of the impulse stages on the horizon. 
  ///
  int N_impulse() const;

  ///
  /// @return Number of the lift stages on the horizon. 
  ///
  int N_lift() const; 

  ///
  /// @return Ideal number of the discretization grids on the horizon. 
  ///
  int N_ideal() const;

  ///
  /// @brief Returns the contact phase of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Contact phase of the time stage. 
  ///
  int contactPhase(const int time_stage) const;

  ///
  /// @brief Returns the contact phase after the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Contact phase after the impulse. 
  ///
  int contactPhaseAfterImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the contact phase after the lift. 
  /// @param[in] lift_index Index of the lift of interest. 
  /// @return Contact phase after the lift. 
  ///
  int contactPhaseAfterLift(const int lift_index) const;

  ///
  /// @brief Returns the impulse index after the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Impulse index after the time stage. 
  ///
  int impulseIndexAfterTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the lift index after the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Lift index after the time stage. 
  ///
  int liftIndexAfterTimeStage(const int time_stage) const;

  ///
  /// @brief Returns the time stage before the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Time stage before the impulse. 
  ///
  int timeStageBeforeImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the time stage after the impulse. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return Time stage after the impulse. 
  ///
  int timeStageAfterImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the time stage before the lift. 
  /// @param[in] lift_index Index of the lift of interest. 
  /// @return Time stage before the lift. 
  ///
  int timeStageBeforeLift(const int lift_index) const;

  ///
  /// @brief Returns the time stage after the lift. 
  /// @param[in] lift_index Index of the lift of interest. 
  /// @return Time stage after the lift. 
  ///
  int timeStageAfterLift(const int lift_index) const;

  ///
  /// @brief Checks wheather the time stage is just before the impulse or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just before the impulse. false if not.
  ///
  bool isTimeStageBeforeImpulse(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just after the impulse or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just after the impulse. false if not.
  ///
  bool isTimeStageAfterImpulse(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just before the lift or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just before the lift. false if not.
  ///
  bool isTimeStageBeforeLift(const int time_stage) const;

  ///
  /// @brief Checks wheather the time stage is just after the lift or not. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return true if the time stage is just after the lift. false if not.
  ///
  bool isTimeStageAfterLift(const int time_stage) const;

  ///
  /// @brief Returns the time of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Time of the time stage of interest.
  ///
  double t(const int time_stage) const;

  ///
  /// @brief Returns the time of the impulse. 
  /// @param[in] impulse_index Index of impulse of interest. 
  /// @return Time of the impulse of interest.
  ///
  double t_impulse(const int impulse_index) const;

  ///
  /// @brief Returns the time of the lift. 
  /// @param[in] lift_index Index of lift of interest. 
  /// @return Time of the lift of interest.
  ///
  double t_lift(const int lift_index) const;

  ///
  /// @brief Returns the time step of the time stage. 
  /// @param[in] time_stage Time stage of interest. 
  /// @return Time step of the time stage of interest.
  ///
  double dt(const int time_stage) const;

  ///
  /// @brief Returns the time step of the auxiliary stage of the impulse. 
  /// @param[in] impulse_index Index of impulse of interest. 
  /// @return Time step of the auxiliary stage of the impulse.
  ///
  double dt_aux(const int impulse_index) const;

  ///
  /// @brief Returns the time step of the auxiliary stage of the lift. 
  /// @param[in] lift_index Index of lift of interest. 
  /// @return Time step of the auxiliary stage of the lift.
  ///
  double dt_lift(const int lift_index) const;

  ///
  /// @brief Returns the ideal time step. 
  /// @return The ideal time step.
  ///
  double dt_ideal() const;

  ///
  /// @brief Checks wheather the STO is enabled for the specified impulse event. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return true if the STO is enabled. false if not.
  ///
  bool isSTOEnabledImpulse(const int impulse_index) const;

  ///
  /// @brief Checks wheather the STO is enabled for the specified lift event. 
  /// @param[in] lift_index Index of the lift of interest. 
  /// @return true if the STO is enabled. false if not.
  ///
  bool isSTOEnabledLift(const int lift_index) const;

  ///
  /// @brief Returns the event index of the specified impulse event. 
  /// @param[in] impulse_index Index of the impulse of interest. 
  /// @return The event index of the specified impulse event.
  ///
  int eventIndexImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the event index of the specified lift event. 
  /// @param[in] lift_index Index of the lift of interest. 
  /// @return The event index of the specified lift event.
  ///
  int eventIndexLift(const int lift_index) const;

  ///
  /// @brief Returns the event type of this discrete event.
  /// @return Event type of this discrete event. 
  ///
  DiscreteEventType eventType(const int event_index) const;

  ///
  /// @brief Checks wheather the optimal control problem is tractable. 
  /// @return true if the optimal control problem is consistent. false if not.
  ///
  bool isFormulationTractable() const;

  ///
  /// @brief Checks wheather the switching times are consistent. 
  /// @return true if the switching times are consistent. false if not.
  ///
  bool isSwitchingTimeConsistent() const;

  ///
  /// @brief Displays the hybrid OCP discretization onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const HybridOCPDiscretization& discretization);

  ///
  /// @brief Minimum step size of the discretization grid. 
  ///
  static constexpr double k_min_dt 
      = std::sqrt(std::numeric_limits<double>::epsilon());

private:
  double T_, dt_ideal_, max_dt_;
  int N_, N_ideal_, N_impulse_, N_lift_, max_events_;
  std::vector<int> contact_phase_index_from_time_stage_, 
                   impulse_index_after_time_stage_, 
                   lift_index_after_time_stage_, time_stage_before_impulse_, 
                   time_stage_before_lift_;
  std::vector<bool> is_time_stage_before_impulse_, is_time_stage_before_lift_,
                    sto_impulse_, sto_lift_;
  std::vector<double> t_, t_impulse_, t_lift_, dt_, dt_aux_, dt_lift_;
  std::vector<DiscreteEventType> event_types_;

  void countDiscreteEvents(
      const std::shared_ptr<ContactSequence>& contact_sequence, const double t);

  void countTimeSteps(const double t);

  void countTimeStages();

  void countContactPhase();

};

} // namespace robotoc

#include "robotoc/hybrid/hybrid_ocp_discretization.hxx"

#endif // ROBOTOC_HYBRID_OCP_DISCRETIZATION_HPP_ 