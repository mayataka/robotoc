#ifndef ROBOTOC_CONTACT_PLANNER_BASE_HPP_
#define ROBOTOC_CONTACT_PLANNER_BASE_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {

#define DEFINE_DEFAULT_CLONE_CONTACT_PLANNER(Derived) \
  std::shared_ptr<ContactPlannerBase> clone() const override { return std::make_shared<Derived>(*this); } 

///
/// @class ContactPlannerBase
/// @brief Base interface of contact planners.
///
class ContactPlannerBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  ContactPlannerBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~ContactPlannerBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ContactPlannerBase(const ContactPlannerBase&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ContactPlannerBase& operator=(const ContactPlannerBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactPlannerBase(ContactPlannerBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactPlannerBase& operator=(ContactPlannerBase&&) noexcept = default;

  ///
  /// @brief Clones this to a shared ptr. 
  ///
  virtual std::shared_ptr<ContactPlannerBase> clone() const = 0;

  ///
  /// @brief Initializes the planner. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  ///
  virtual void init(const Eigen::VectorXd& q) = 0;

  ///
  /// @brief Plans the foot steps in the MPC (receding-horizon) fashion. 
  /// @param[in] t Initial time.
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] contact_status Initial contact status.
  /// @param[in] planning_steps Number of planning steps. Must be non-negative.
  /// @return True if the planning is succeeded. False if not.
  /// @remark The linear and angular velocities of the floating base are assumed
  /// to be expressed in the body local coordinate.
  /// @remark The implementation must follow: step=0: previous step, 
  /// step=1: initial step (specified as q and contact_status).
  ///
  virtual bool plan(const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, 
                    const ContactStatus& contact_status,
                    const int planning_steps) = 0;

  ///
  /// @brief Gets the contact placements of a specified step. 
  /// @param[in] step Step of interest.
  /// @return const reference to the contact placements of a specified step. 
  /// @remark step=0: previous step, step=1: initial step.
  ///
  virtual const aligned_vector<SE3>& contactPlacements(const int step) const = 0;

  ///
  /// @brief Gets the contact positions of a specified step. 
  /// @param[in] step Step of interest.
  /// @return const reference to the contact positions of a specified step. 
  /// @remark step=0: previous step, step=1: initial step.
  ///
  virtual const std::vector<Eigen::Vector3d>& contactPositions(const int step) const = 0;

  ///
  /// @brief Gets the contact surfaces of a specified step. 
  /// @param[in] step Step of interest.
  /// @return const reference to the contact surfaces of a specified step. 
  /// @remark step=0: previous step, step=1: initial step.
  ///
  virtual const std::vector<Eigen::Matrix3d>& contactSurfaces(const int step) const = 0;

  ///
  /// @brief Gets the CoM position of a specified step. 
  /// @param[in] step Step of interest.
  /// @return const reference to CoM position of a specified step. 
  /// @remark step=0: previous step, step=1: initial step.
  ///
  virtual const Eigen::Vector3d& CoM(const int step) const = 0;

  ///
  /// @brief Gets the rotation matrix of the base at a specified step. 
  /// @param[in] step Step of interest.
  /// @return const reference to rotation matrix of the base at a specified step. 
  /// @remark step=0: previous step, step=1: initial step.
  ///
  virtual const Eigen::Matrix3d& R(const int step) const = 0;

  ///
  /// @brief Gets the size of each planned quantities. 
  /// @return Size of planning.
  ///
  virtual int size() const = 0;

  void disp(std::ostream& os) const;

private:

};

} // namespace robotoc 

#endif // ROBOTOC_CONTACT_PLANNER_BASE_HPP_ 