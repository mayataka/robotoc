#ifndef IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HPP_
#define IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HPP_ 

#include <vector>

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/hybrid/discrete_event.hpp"


namespace idocp {

///
/// @class ContactSequencePrimitive 
/// @brief Primiteive implementation of the contact sequence, i.e., sequence of  
/// contact status over the horizon. 
///
class ContactSequencePrimitive {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  ///
  ContactSequencePrimitive(const Robot& robot, const int N);

  ///
  /// @brief Default constructor. 
  ///
  ContactSequencePrimitive();

  ///
  /// @brief Destructor. 
  ///
  ~ContactSequencePrimitive();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactSequencePrimitive(const ContactSequencePrimitive&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ContactSequencePrimitive& operator=(const ContactSequencePrimitive&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactSequencePrimitive(ContactSequencePrimitive&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactSequencePrimitive& operator=(ContactSequencePrimitive&&) noexcept = default;

  ///
  /// @brief Set the contact status over all of the time stages uniformly. Also, 
  /// disable discrete events over all of the time stages.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatusUniformly(const ContactStatus& contact_status);

  ///
  /// @brief Set the discrete event. Contact status after discrete event is also
  /// uniformly changed by discrete_event. 
  /// @param[in] discrete_event Discrete event.
  /// @param[in] event_time_stage Event time stage.
  ///
  void setDiscreteEvent(const DiscreteEvent& discrete_event,
                        const int event_time_stage);

  ///
  /// @brief Shift the discrete event. 
  /// @param[in] event_time_stage Even time stage where the discrete event exists.
  /// ContactSequencePrimitive::existDiscreteEvent() must be true.
  /// @param[in] shifted_event_time_stage Event time stage where the discrete 
  /// event shifts.
  ///
  void shiftDiscreteEvent(const int event_time_stage, 
                          const int shifted_event_time_stage);

  ///
  /// @brief Shift the discrete event beyond initial.
  /// @param[in] event_time_stage Event time stage where the discrete event exists.
  /// ContactSequencePrimitive::existDiscreteEvent() must be true.
  ///
  void shiftDiscreteEventBeyondInitial(const int event_time_stage);

  ///
  /// @brief Shift the discrete event beyond terminal.
  /// @param[in] event_time_stage Event time stage where the discrete event exists.
  /// ContactSequencePrimitive::existDiscreteEvent() must be true.
  ///
  void shiftDiscreteEventBeyondTerminal(const int event_time_stage);

  ///
  /// @brief Getter of the contact status. 
  /// @param[in] time_stage Time stage of interest.
  /// @return const reference to the contact status.
  ///
  const ContactStatus& contactStatus(const int time_stage) const;

  ///
  /// @brief Getter of the impulse status. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// ContactSequencePrimitive::existImpulse() must be true.
  /// @return const reference to the impulse status.
  ///
  const ImpulseStatus& impulseStatus(const int event_time_stage) const;

  ///
  /// @brief Getter of the event time. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// ContactSequencePrimitive::existDiscreteEvent() must be true.
  /// @return Event time.
  ///
  double eventTime(const int event_time_stage) const;

  ///
  /// @brief Setter of the event time. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// ContactSequencePrimitive::existDiscreteEvent() must be true.
  /// @param[in] event_time Event time.
  ///
  void setEventTime(const int event_time_stage, const double event_time);

  ///
  /// @brief Check wheather a discrete event exists at a time stage. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// @return True if a discrete event exists at time_stage. False otherwise.
  ///
  bool existDiscreteEvent(const int event_time_stage) const;

  ///
  /// @brief Check wheather a impulse exists at a time stage. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// @return True if a impulse exists at time_stage. False otherwise.
  ///
  bool existImpulse(const int event_time_stage) const;

  ///
  /// @brief Check wheather a lift exists at a time stage. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// @return True if a lift exists at time_stage. False otherwise.
  ///
  bool existLift(const int event_time_stage) const;

  ///
  /// @brief Check wheather a lift exists and impulse does not at a time stage. 
  /// @param[in] event_time_stage Event time stage of interest.
  /// @return True if a lift exists and impulse does not at time_stage. 
  /// False otherwise.
  ///
  bool existOnlyLift(const int event_time_stage) const;

private:
  int N_;
  std::vector<ContactStatus> contact_sequence_;
  std::vector<DiscreteEvent> discrete_event_sequence_;

};

} // namespace idocp 

#include "idocp/hybrid/contact_sequence_primitive.hxx"

#endif // IDOCP_CONTACT_SEQUENCE_PRIMITIVE_HPP_ 