#ifndef IDOCP_DISCRETE_EVENT_HPP_
#define IDOCP_DISCRETE_EVENT_HPP_

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class DiscreteEvent
/// @brief Discrete event composed by impulse and lift.
///
class DiscreteEvent {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_point_contacts Maximum number of the point contacts.
  ///
  DiscreteEvent(const int max_point_contacts);

  ///
  /// @brief Default constructor. 
  ///
  DiscreteEvent();

  ///
  /// @brief Destructor. 
  ///
  ~DiscreteEvent();

  ///
  /// @brief Default copy constructor. 
  ///
  DiscreteEvent(const DiscreteEvent&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  DiscreteEvent& operator=(const DiscreteEvent&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  DiscreteEvent(DiscreteEvent&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  DiscreteEvent& operator=(DiscreteEvent&&) noexcept = default;

  ///
  /// @brief Returns const reference to impulse status. 
  /// @return const reference to impulse status. 
  ///
  const ImpulseStatus& impulseStatus() const;

  ///
  /// @brief Returns true if this object has discrete event. Returns false if not.
  /// @return true if this object has discrete event. false if not.
  ///
  bool hasDiscreteEvent() const;

  ///
  /// @brief Returns true if this object has impulse. Returns false if not.
  /// @return true if this object has impulse. false if not.
  ///
  bool hasImpulse() const;

  ///
  /// @brief Returns true if this object has lift. Returns false if not.
  /// @return true if this object has lift. false if not.
  ///
  bool hasLift() const;

  ///
  /// @brief Returns the time of the discrete event occurs.
  /// @return The time of the discrete event occurs.
  ///
  double eventTime() const;

  ///
  /// @brief Act on the contact status and change it according to this discrete 
  /// event.
  /// @param[in] Contact status. Must be consistent with this discrete event. 
  ///
  void act(ContactStatus& contact_status) const;

  ///
  /// @brief Set the contact status from two sequential contact status.
  /// @param[in] contact_status_before Contact status before this discrete event. 
  /// @param[in] contact_status_after Contact status after this discrete event. 
  ///
  void setDiscreteEvent(const ContactStatus& contact_status_before, 
                        const ContactStatus& contact_status_after);

  ///
  /// @brief set the time of the discrete event.
  ///
  void setEventTime(const double time);

  ///
  /// @brief Disable discrete event.
  ///
  void disableDiscreteEvent();

  ///
  /// @brief Return the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int max_point_contacts() const;

private:
  ContactStatus contact_status_before_, contact_status_after_;
  ImpulseStatus impulse_status_;
  int max_point_contacts_;
  bool has_impulse_, has_lift_;
  double time_;

};

} // namespace idocp

#include "idocp/ocp/discrete_event.hxx"

#endif // IDOCP_DISCRETE_EVENT_HPP_