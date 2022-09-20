#ifndef ROBOTOC_ROBOT_MODEL_INFO_HPP_
#define ROBOTOC_ROBOT_MODEL_INFO_HPP_

#include <string>
#include <vector>

#include "robotoc/robot/contact_model_info.hpp"

namespace robotoc {

///
/// @enum BaseJointType
/// @brief Types of the base joints of robots
///
enum class BaseJointType {
  FixedBase,
  FloatingBase
};

///
/// @class RobotModelInfo
/// @brief Info of a robot model. 
///
struct RobotModelInfo {
  ///
  /// @brief Construct a robot model info.
  /// @param[in] urdf_path Path to the URDF file.
  /// @param[in] base_joint_type Type of the base joint. 
  /// @param[in] point_contacts Info of point contacts.
  /// @param[in] surface_contacts Info of surface contacts.
  /// @param[in] contact_inv_damping Damping paramter in matrix inversion of the 
  /// contact-consistent forward dynamics. Default is 0.0.
  ///
  RobotModelInfo(const std::string& urdf_path, 
                 const BaseJointType base_joint_type,
                 const std::vector<ContactModelInfo>& point_contacts,
                 const std::vector<ContactModelInfo>& surface_contacts,
                 const double contact_inv_damping=0.0);

  ///
  /// @brief Default constructor. 
  ///
  RobotModelInfo() = default;

  ///
  /// @brief Default destructor. 
  ///
  ~RobotModelInfo() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  RobotModelInfo(const RobotModelInfo&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  RobotModelInfo& operator=(const RobotModelInfo&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RobotModelInfo(RobotModelInfo&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RobotModelInfo& operator=(RobotModelInfo&&) noexcept = default;

  ///
  /// @brief Path to the URDF file.
  ///
  std::string urdf_path;

  ///
  /// @brief Type of the base joint. Default is BaseJointType::FixedBase.
  ///
  BaseJointType base_joint_type = BaseJointType::FixedBase;

  ///
  /// @brief Info of point contacts.
  ///
  std::vector<ContactModelInfo> point_contacts;

  ///
  /// @brief Info of surface contacts.
  ///
  std::vector<ContactModelInfo> surface_contacts;

  ///
  /// @brief Damping paramter in matrix inversion of the contact-consistent 
  /// forward dynamics. 1e-12 works well for two surface contacts. 
  /// Must be non-negative. Default is 0.
  ///
  double contact_inv_damping = 0.0;

  ///
  /// @brief Creates a simple robot manipulator model info.
  /// @param[in] urdf_path Path to the URDF file.
  ///
  static RobotModelInfo Manipulator(const std::string& urdf_path);

  ///
  /// @brief Creates a simple quadrped robot model info.
  /// @param[in] urdf_path Path to the URDF file.
  /// @param[in] point_contacts Info of point contacts. 
  ///
  static RobotModelInfo Quadruped(
      const std::string& urdf_path, 
      const std::vector<ContactModelInfo>& point_contacts);

  ///
  /// @brief Creates a simple humanoid robot model info.
  /// @param[in] urdf_path Path to the URDF file.
  /// @param[in] surface_contacts Info of surface contacts.
  ///
  static RobotModelInfo Humanoid(
      const std::string& urdf_path, 
      const std::vector<ContactModelInfo>& surface_contacts);

};

} // namespace robotoc

#endif // ROBOTOC_ROBOT_MODEL_INFO_HPP_