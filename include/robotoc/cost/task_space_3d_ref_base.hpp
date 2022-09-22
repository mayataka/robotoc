#ifndef ROBOTOC_TASK_SPACE_3D_REF_BASE_HPP_
#define ROBOTOC_TASK_SPACE_3D_REF_BASE_HPP_

#include "Eigen/Core"

#include "robotoc/ocp/grid_info.hpp"


namespace robotoc {

///
/// @class TaskSpace3DRefBase
/// @brief Base class of reference task space position. 
///
class TaskSpace3DRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TaskSpace3DRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TaskSpace3DRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TaskSpace3DRefBase(const TaskSpace3DRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TaskSpace3DRefBase& operator=(const TaskSpace3DRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TaskSpace3DRefBase(TaskSpace3DRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TaskSpace3DRefBase& operator=(TaskSpace3DRefBase&&) noexcept = default;

  ///
  /// @brief Computes the reference task-space position. 
  /// @param[in] grid_info Grid info.
  /// @param[in] ref_3d Reference position. Size is 3.
  ///
  virtual void updateRef(const GridInfo& grid_info, 
                         Eigen::VectorXd& ref_3d) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] grid_info Grid info.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const GridInfo& grid_info) const = 0;
};

} // namespace robotoc


#endif // ROBOTOC_TASK_SPACE_3D_REF_BASE_HPP_