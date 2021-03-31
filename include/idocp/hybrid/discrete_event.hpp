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
  /// @brief Set the contact status from two sequential contact status.
  /// @param[in] pre_contact_status Contact status before this discrete event. 
  /// @param[in] post_contact_status Contact status after this discrete event. 
  ///
  DiscreteEvent(const ContactStatus& pre_contact_status, 
                const ContactStatus& post_contact_status);

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
  /// @brief Returns const reference to impulse status. 
  /// @return const reference to impulse status. 
  ///
  const ImpulseStatus& impulseStatus() const;

  ///
  /// @brief Gets the contact status before this discrete event.
  /// @return const reference to the pre contact status.
  ///
  const ContactStatus& preContactStatus() const;

  ///
  /// @brief Gets the contact status after this discrete event.
  /// @return const reference to the post contact status.
  ///
  const ContactStatus& postContactStatus() const;

  ///
  /// @brief Sets the contact status from two sequential contact status.
  /// @param[in] pre_contact_status Contact status before this discrete event. 
  /// @param[in] post_contact_status Contact status after this discrete event. 
  ///
  void setDiscreteEvent(const ContactStatus& pre_contact_status, 
                        const ContactStatus& post_contact_status);

  ///
  /// @brief Sets a contact point.
  /// @param[in] contact_index Index of the contact.
  /// @param[in] contact_point Contact point.
  ///
  void setContactPoint(const int contact_index, 
                       const Eigen::Vector3d& contact_point);

  ///
  /// @brief Sets contact points.
  /// @param[in] contact_points Contact points. Size must be 
  /// ImpulseStatus::maxPointContacts().
  ///
  void setContactPoints(const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Returns the maximum number of the contacts.
  /// @return The maximum number of the contacts. 
  ///
  int maxPointContacts() const;

private:
  ContactStatus pre_contact_status_, post_contact_status_;
  ImpulseStatus impulse_status_;
  int max_point_contacts_;
  bool exist_impulse_, exist_lift_;

};

} // namespace idocp

#include "idocp/hybrid/discrete_event.hxx"

#endif // IDOCP_DISCRETE_EVENT_HPP_