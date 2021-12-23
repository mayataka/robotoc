#ifndef ROBOTOC_CONTACT_FRAMES_HPP_
#define ROBOTOC_CONTACT_FRAMES_HPP_

#include <vector>


namespace robotoc {

///
/// @class ContactFrames
/// @brief Collection of the frames involving the point and surface contacts.
///
class ContactFrames {
public:
  ///
  /// @brief Constructor.
  /// @param[in] point_contact_frames Point contact frames.
  /// @param[in] surface_contact_frames Surface contact frames.
  ///
  ContactFrames(const std::vector<int>& point_contact_frames, 
                const std::vector<int>& surface_contact_frames);

  ///
  /// @brief Default constructor. 
  ///
  ContactFrames();

  ///
  /// @brief Destructor. 
  ///
  ~ContactFrames();

  ///
  /// @brief Use default copy constructor. 
  ///
  ContactFrames(const ContactFrames&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  ContactFrames& operator=(const ContactFrames&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  ContactFrames(ContactFrames&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  ContactFrames& operator=(ContactFrames&&) noexcept = default;

  /// 
  /// @brief Collection of the frames involving the point contacts.
  /// 
  std::vector<int> point_contact_frames;

  /// 
  /// @brief Collection of the frames involving the surface contacts.
  /// 
  std::vector<int> surface_contact_frames;
};

} // namespace robotoc

#endif // ROBOTOC_CONTACT_FRAMES_HPP_ 