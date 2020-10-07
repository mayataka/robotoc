#ifndef IDOCP_CONTACT_SEQUENCE_HPP_
#define IDOCP_CONTACT_SEQUENCE_HPP_

#include <vector>

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


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
  ///
  ContactSequence(const Robot& robot, const int N);

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
  /// @brief Activate a contact over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] time_stage_begin Start time stage. 
  /// @param[in] time_stage_end Last time stage. 
  ///
  void activateContact(const int contact_index, const int time_stage_begin, 
                       const int time_stage_end);

  ///
  /// @brief Deactivate a contact over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @param[in] time_stage_begin Start time stage. 
  /// @param[in] time_stage_end Last time stage. 
  ///
  void deactivateContact(const int contact_index, const int time_stage_begin, 
                         const int time_stage_end);

  ///
  /// @brief Activate contacts over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] time_stage_begin Start time stage. 
  /// @param[in] time_stage_end Last time stage. 
  ///
  void activateContacts(const std::vector<int>& contact_indices, 
                        const int time_stage_begin, const int time_stage_end);

  ///
  /// @brief Deactivate contacts over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] contact_indices Indices of contacts of interedted. 
  /// @param[in] time_stage_begin Start time stage. 
  /// @param[in] time_stage_end Last time stage. 
  ///
  void deactivateContacts(const std::vector<int>& contact_indices, 
                          const int time_stage_begin, 
                          const int time_stage_end);

  ///
  /// @brief Deactivate contacts over specified time steps 
  /// (from time_stage_begin to time_stage_end). 
  /// @param[in] time_stage Last time stage. 
  ///
  const ContactStatus& contactStatus(const int time_stage) const;

private:
  int max_point_contacts_, N_;
  std::vector<ContactStatus> contact_sequence_;

};

} // namespace idocp 

#include "idocp/ocp/contact_sequence.hxx"

#endif // IDOCP_CONTACT_SEQUENCE_HPP_