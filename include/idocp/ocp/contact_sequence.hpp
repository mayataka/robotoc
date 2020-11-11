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
  /// horizon. 
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
  /// @brief Set all of the contact status uniformly. Also, disable all of the 
  /// discrete events.
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
  /// @brief Shift the discrete event. 
  /// @param[in] time_stage Time stage where the discrete event of interest 
  /// occurs.
  /// @param[in] event_time Event time after shifted.
  ///
  void shiftDiscreteEvent(const int time_stage, const double event_time);

  ///
  /// @brief Shift the discrete event. 
  /// @param[in] time_stage Time stage where the discrete event of interest 
  /// occurs.
  ///
  void shiftDiscreteEventBeyondInitial(const int time_stage);

  ///
  /// @brief Shift the discrete event. 
  /// @param[in] time_stage Time stage where the discrete event of interest 
  /// occurs.
  ///
  void shiftDiscreteEventBeyondTerminal(const int time_stage);

  ///
  /// @brief Getter of the contact status. 
  /// @param[in] time_stage Time stage of interest.
  /// @return const reference to the contact status.
  ///
  const ContactStatus& contactStatus(const int time_stage) const;

  ///
  /// @brief Getter of the impulse status. 
  /// @param[in] time_stage Time stage of interest.
  /// @return const reference to the impulse status.
  ///
  const ImpulseStatus& impulseStatus(const int time_stage) const;

  ///
  /// @brief Getter of the event time. 
  /// @param[in] time_stage Time stage of interest.
  /// @return Event time.
  ///
  double eventTime(const int time_stage) const;

  ///
  /// @brief Returns true if discrete event exists over time_stage. Returns 
  /// false if not.
  /// @param[in] time_stage Time stage of interested.
  /// @return true if discrete event exists over time_stage. false if not.
  ///
  bool existDiscreteEvent(const int time_stage) const;

  ///
  /// @brief Returns true if impulse exists over time_stage. Returns false if
  /// not.
  /// @param[in] time_stage Time stage of interested.
  /// @return true if impulse exists over time_stage. false if not.
  ///
  bool existImpulse(const int time_stage) const;

  ///
  /// @brief Returns the index of the impulse at time_stage. Returned value is
  /// valid only if ContactSequence::existImpulse(time_stage) is true. 
  /// Otherwise return -1.
  /// @param[in] time_stage Time stage of interested.
  /// @return Index of impulse.
  ///
  int impulseIndex(const int time_stage) const;

  ///
  /// @brief Returns true if lift exists over time_stage. Returns false if not.
  /// @param[in] time_stage Time stage of interested.
  /// @return true if lift exists over time_stage. false if not.
  ///
  bool existLift(const int time_stage) const;

  ///
  /// @brief Returns the index of the lift at time_stage. Returned value is 
  /// valid only if ContactSequence::existLift(time_stage) is true.
  /// Otherwise return -1.
  /// @param[in] time_stage Time stage of interested.
  /// @return Index of lift.
  ///
  int liftIndex(const int time_stage) const;

  ///
  /// @brief Returns the number of impulse.
  /// @return Number of impulse.
  ///
  int numImpulse() const;

  ///
  /// @brief Returns the number of lift.
  /// @return Number of lift.
  ///
  int numLift() const;

  ///
  /// @brief Returns the time stage of the impulse having impulse_index.  
  /// Returned value is valid only if impulse_index is less than 
  /// ContactSequence::numImpoulse(). Otherwise return -1.
  /// @param[in] impulse_index Index of the impulse of interested.
  /// @return Time stage of impulse.
  ///
  int impulseStage(const int impulse_index) const;

  ///
  /// @brief Returns the time stage of the lift having lift_index.  
  /// Returned value is valid only if lift_index is less than 
  /// ContactSequence::numLift(). Otherwise return -1.
  /// @param[in] lift_index Index of the impulse of interested.
  /// @return Time stage of lift.
  ///
  int liftStage(const int lift_index) const;

private:
  int N_;
  double T_, dtau_;
  std::vector<ContactStatus> contact_sequence_;
  std::vector<DiscreteEvent> discrete_event_sequence_;
  std::vector<int> impulse_index_, lift_index_, impulse_stage_, lift_stage_;

  int timeStageFromTime(const double time) const;

  void countAll();

  void countImpulse();

  void countLift();

  void countImpulseStage();

  void countLiftStage();
};

} // namespace idocp 

#include "idocp/ocp/contact_sequence.hxx"

#endif // IDOCP_CONTACT_SEQUENCE_HPP_