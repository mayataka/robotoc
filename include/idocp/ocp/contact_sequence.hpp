#ifndef IDOCP_CONTACT_SEQUENCE_HPP_
#define IDOCP_CONTACT_SEQUENCE_HPP_

#include <vector>

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/discrete_event.hpp"


namespace idocp {

///
/// @class ContactSequence 
/// @brief Contact sequence, i.e., sequence of contact status over the 
/// horizon. Provides the formulation of the optimal control problem 
/// with impulse and lift.
///
class ContactSequence {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  ///
  ContactSequence(const Robot& robot, const double T, const int N);

  ///
  /// @brief Default constructor. 
  ///
  ContactSequence();

  ///
  /// @brief Destructor. 
  ///
  ~ContactSequence();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactSequence(const ContactSequence&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ContactSequence& operator=(const ContactSequence&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactSequence(ContactSequence&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactSequence& operator=(ContactSequence&&) noexcept = default;

  ///
  /// @brief Set all of the contact status uniformly. Also, disable all of 
  /// discrete events, i.e., impulse and lift.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatusUniformly(const ContactStatus& contact_status);

  ///
  /// @brief Set the discrete event. Contact status after discrete event is also
  /// uniformly changed by discrete_event. To determine the event time and 
  /// stage, DiscreteEvent::eventTime() is internally called. Hence, set the 
  /// time by DiscreteEvent::setEventTime() before calling this function.
  /// @param[in] discrete_event Discrete event.
  ///
  void setDiscreteEvent(const DiscreteEvent& discrete_event);

  ///
  /// @brief Shift the impulse. 
  /// @param[in] impulse_index Impulse index. Must be less than 
  /// ContactSequence::totalNumImpulseStages().
  /// @param[in] impulse_time Impulse time after shifted.
  ///
  void shiftImpulse(const int impulse_index, const double impulse_time);

  ///
  /// @brief Shift the discrete event. 
  /// @param[in] lift_index Lift index. Must be less than 
  /// ContactSequence::totalNumLiftStages().
  /// @param[in] lift_time Lift time after shifted.
  ///
  void shiftLift(const int lift_index, const double lift_time);

  ///
  /// @brief Getter of the contact status. 
  /// @param[in] time_stage Time stage of interest. The discrete event of  
  /// interest occurs between time_stage and time_stage+1.
  /// @return const reference to the contact status.
  ///
  const ContactStatus& contactStatus(const int time_stage) const;

  ///
  /// @brief Getter of the impulse status. 
  /// ContactSequence::totalNumImpulseStages() must be positive.
  /// @param[in] impulse_index Impulse index. Must be less than 
  /// ContactSequence::totalNumImpulseStages().
  /// @return const reference to the impulse status.
  ///
  const ImpulseStatus& impulseStatus(const int impulse_index) const;

  ///
  /// @brief Getter of the impulse time. 
  /// @param[in] impulse_index Impulse index. Must be less than 
  /// ContactSequence::totalNumImpulseStages().
  /// @return Impulse time.
  ///
  double impulseTime(const int impulse_index) const;

  ///
  /// @brief Getter of the lift time. 
  /// @param[in] lift_index Lift index. Must be less than 
  /// ContactSequence::totalNumLiftStages().
  /// @return Lift time.
  ///
  double liftTime(const int lift_index) const;

  ///
  /// @brief Returns the number of impulse stages before the time_stage. 
  /// @param[in] time_stage Time stage of interested.
  /// @return Number of impulse stages.
  ///
  int numImpulseStages(const int time_stage) const;

  ///
  /// @brief Returns the number of lift stages before the time_stage. 
  /// @param[in] time_stage Time stage of interested.
  /// @return Number of lift stages.
  ///
  int numLiftStages(const int time_stage) const;

  ///
  /// @brief Returns the total number of impulse stages over the horizon.
  /// @return Total number of impulse stages over the horizon.
  ///
  int totalNumImpulseStages() const;

  ///
  /// @brief Returns the total number of lift stages over the horizon.
  /// @return Total number of lift stages over the horizon.
  ///
  int totalNumLiftStages() const;

  ///
  /// @brief Returns the time stage just before impulse having impulse_index.  
  /// @param[in] impulse_index Impulse index. Must be less than 
  /// ContactSequence::totalNumImpulseStages().
  /// @return Time stage before the impulse.
  ///
  int timeStageBeforeImpulse(const int impulse_index) const;

  ///
  /// @brief Returns the time stage just before lift having lift_index.  
  /// @param[in] lift_index Lift index. Must be less than 
  /// ContactSequence::totalNumLiftStages().
  /// @return Time stage before the lift.
  ///
  int timeStageBeforeLift(const int lift_index) const;

private:
  int N_;
  double T_, dtau_;
  std::vector<ContactStatus> contact_sequence_;
  std::vector<DiscreteEvent> discrete_event_sequence_;
  std::vector<int> num_impulse_stages_, num_lift_stages_,  
                   impulse_stage_, lift_stage_;

  void shiftDiscreteEvent(const int event_index, const double event_time);

  void shiftDiscreteEventBeyondInitial(const int time_stage);

  void shiftDiscreteEventBeyondTerminal(const int time_stage);

  bool existDiscreteEvent(const int time_stage) const;

  bool existImpulse(const int time_stage) const;

  bool existLift(const int time_stage) const;

  int timeStageFromContinuousTime(const double time) const;

  void countAll();

  void countNumImpulseStaeges();

  void countNumLiftStages();

  void countImpulseStage();

  void countLiftStage();
};

} // namespace idocp 

#include "idocp/ocp/contact_sequence.hxx"

#endif // IDOCP_CONTACT_SEQUENCE_HPP_