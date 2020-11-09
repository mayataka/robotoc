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
  ///
  void setContactStatusUniformly(const ContactStatus& contact_status);

  void setDiscreteEvent(const DiscreteEvent& discrete_event);

  void setDiscreteEvent(const int time_stage, 
                        const DiscreteEvent& discrete_event);

  void shiftDiscreteEvent(const int time_stage, const double shift_event_time);

  void shiftDiscreteEvent(const int time_stage, const int shit_time_stages);

  void removeDiscreteEvent(const int time_stage);

  ///
  /// @brief Deactivate contacts over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] time_stage Last time stage. 
  ///
  const ContactStatus& contactStatus(const int time_stage) const;

  const ImpulseStatus& impulseStatus(const int time_stage) const;

  bool hasDiscreteEvent(const int time_stage) const;

  bool hasImpulse(const int time_stage) const;

  bool hasLift(const int time_stage) const;

private:
  int max_point_contacts_, N_;
  double T_, dtau_;
  std::vector<ContactStatus> contact_sequence_;
  std::vector<DiscreteEvent> discrete_event_sequence_;

  void setContactSequenceFromDiscreteEvent(const int time_stage_begin);
};

} // namespace idocp 

#include "idocp/ocp/contact_sequence.hxx"

#endif // IDOCP_CONTACT_SEQUENCE_HPP_