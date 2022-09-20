#ifndef ROBOTOC_CONTACT_MODEL_INFO_HPP_
#define ROBOTOC_CONTACT_MODEL_INFO_HPP_

#include <string>

namespace robotoc {

///
/// @class ContactModelInfo
/// @brief Info of a contact model. 
///
struct ContactModelInfo {
  ///
  /// @brief Construct a contact model info.
  /// @param[in] frame Name of the contact frame.
  /// @param[in] baumgarte_time_step Time step parameter of the Baumgarte's 
  /// stabilization method. Must be positive.
  ///
  ContactModelInfo(const std::string& frame, const double baumgarte_time_step);

  ///
  /// @brief Construct a contact model info.
  /// @param[in] frame Name of the contact frame.
  /// @param[in] baumgarte_position_gain The position gain of the Baumgarte's 
  /// stabilization method. Must be non-negative.
  /// @param[in] baumgarte_velocity_gain The velocity gain of the Baumgarte's 
  /// stabilization method. Must be non-negative.
  ///
  ContactModelInfo(const std::string& frame, 
                   const double baumgarte_position_gain,
                   const double baumgarte_velocity_gain);

  ///
  /// @brief Default constructor. 
  ///
  ContactModelInfo() = default;

  ///
  /// @brief Default destructor. 
  ///
  ~ContactModelInfo() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ContactModelInfo(const ContactModelInfo&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ContactModelInfo& operator=(const ContactModelInfo&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactModelInfo(ContactModelInfo&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactModelInfo& operator=(ContactModelInfo&&) noexcept = default;

  ///
  /// @brief Name of the contact frame.
  ///
  std::string frame;

  ///
  /// @brief The position gain of the Baumgarte's stabilization method. 
  /// Default is 0.0.
  ///
  double baumgarte_position_gain = 0.0;

  ///
  /// @brief The velocity gain of the Baumgarte's stabilization method. 
  /// Default is 0.0.
  ///
  double baumgarte_velocity_gain = 0.0;
};

} // namespace robotoc

#endif // ROBOTOC_CONTACT_MODEL_INFO_HPP_