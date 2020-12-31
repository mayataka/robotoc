#ifndef IDOCP_COLLISION_CHECKER_HPP_
#define IDOCP_COLLISION_CHECKER_HPP_

#include <vector>
#include "Eigen/Core"
#include "idocp/robot/robot.hpp"

namespace idocp {

///
/// @class CollisionChecker
/// @brief A very simple collision checker.
///
class CollisionChecker {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  CollisionChecker(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  CollisionChecker();

  ///
  /// @brief Destructor. 
  ///
  ~CollisionChecker();

  ///
  /// @brief Default copy constructor. 
  ///
  CollisionChecker(const CollisionChecker&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  CollisionChecker& operator=(const CollisionChecker&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CollisionChecker(CollisionChecker&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CollisionChecker& operator=(CollisionChecker&&) noexcept = default;

  ///
  /// @brief Checks collisions. Internally, computes the frame positions for 
  /// each contact frames.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] q Configuration of the robot.
  /// @return If the frames have contacts or not.
  ///
  std::vector<bool> check(Robot& robot, const Eigen::VectorXd& q);

  ///
  /// @brief Get the contact frame positions. Call CollisionChecker::check()
  /// before calling this function.
  /// @return const reference to the contact frame positions.
  ///
  const std::vector<Eigen::Vector3d>& contactFramePositions() const;

  ///
  /// @brief Print the contact frame positions onto the console.
  ///
  void printContactFramePositions() const;

private:
  std::vector<Eigen::Vector3d> contact_points_;
  std::vector<int> contact_frames_indices_;

};

} // namespace idocp

#include "idocp/hybrid/collision_checker.hxx"

#endif // IDOCP_COLLISION_CHECKER_HPP_ 