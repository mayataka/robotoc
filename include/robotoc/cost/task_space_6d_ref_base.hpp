#ifndef ROBOTOC_TASK_SPACE_6D_REF_BASE_HPP_
#define ROBOTOC_TASK_SPACE_6D_REF_BASE_HPP_

#include "Eigen/Core"

#include "robotoc/hybrid/grid_info.hpp"
#include "robotoc/robot/se3.hpp"


namespace robotoc {

///
/// @class TaskSpace6DRefBase
/// @brief Base class of reference task space placement (position and orientation). 
///
class TaskSpace6DRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TaskSpace6DRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TaskSpace6DRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TaskSpace6DRefBase(const TaskSpace6DRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TaskSpace6DRefBase& operator=(const TaskSpace6DRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TaskSpace6DRefBase(TaskSpace6DRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TaskSpace6DRefBase& operator=(TaskSpace6DRefBase&&) noexcept = default;

  ///
  /// @brief Computes the reference task-space placement. 
  /// @param[in] grid_info Grid info.
  /// @param[in] ref_6d Reference placement.
  ///
  virtual void updateRef(const GridInfo& grid_info, SE3& ref_6d) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] grid_info Grid info.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const GridInfo& grid_info) const = 0;
};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_6D_REF_BASE_HPP_