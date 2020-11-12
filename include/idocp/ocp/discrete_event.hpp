#ifndef IDOCP_DISCRETE_EVENT_HPP_
#define IDOCP_DISCRETE_EVENT_HPP_

#include "idocp/robot/robot.hpp"
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
  /// @brief Constructor. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  DiscreteEvent(const Robot& robot);

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
  /// @brief Returns true if this discrete event exists. Returns false if not.
  /// @return true if this discrete event exists. false if not.
  ///
  bool existDiscreteEvent() const;

  ///
  /// @brief Returns true if impulse exists in this discrete event. Returns 
  /// false if not.
  /// @return true if impulse exists. false if not.
  ///
  bool existImpulse() const;

  ///
  /// @brief Returns true if lift exists in this discrete event. Returns false 
  /// if not.
  /// @return true if lift exists. false if not.
  ///
  bool existLift() const;

  ///
  /// @brief Returns the time of the discrete event occurs.
  /// @return The time of the discrete event occurs.
  ///
  double eventTime() const;

  ///
  /// @brief Change the input contact status as is before this discrete event.
  /// @param[in] contact_status Contact status. Must be consistent with this 
  /// discrete event. 
  ///
  void act(ContactStatus& contact_status) const;

  ///
  /// @brief Change the input contact status as is after this discrete event.
  /// @param[in] contact_status Contact status. Must be consistent with this 
  /// discrete event. 
  ///
  void actInv(ContactStatus& contact_status) const;

  ///
  /// @brief Check whether the input contact_status is consistent with the 
  /// pre contact status of this discrete event.
  /// @param[in] contact_status Contact status. 
  /// @return true if input contact_status is consistent with the pre contact 
  /// status of this discrete event. false if not.
  ///
  bool isConsisitentWithPreContactStatus(
      const ContactStatus& contact_status) const;

  ///
  /// @brief Check whether the input contact_status is consistent with the 
  /// post contact status of this discrete event.
  /// @param[in] contact_status Contact status. 
  /// @return true if input contact_status is consistent with the post contact 
  /// status of this discrete event. false if not.
  ///
  bool isConsisitentWithPostContactStatus(
      const ContactStatus& contact_status) const;

  ///
  /// @brief Set the contact status from two sequential contact status.
  /// @param[in] pre_contact_status Contact status before this discrete event. 
  /// @param[in] post_contact_status Contact status after this discrete event. 
  ///
  void setDiscreteEvent(const ContactStatus& pre_contact_status, 
                        const ContactStatus& post_contact_status);

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
  ContactStatus pre_contact_status_, post_contact_status_;
  ImpulseStatus impulse_status_;
  int max_point_contacts_;
  bool exist_impulse_, exist_lift_;
  double time_;

};

} // namespace idocp

#include "idocp/ocp/discrete_event.hxx"

#endif // IDOCP_DISCRETE_EVENT_HPP_