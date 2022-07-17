#ifndef ROBOTOC_CONFIGURATION_SPACE_REF_HPP_
#define ROBOTOC_CONFIGURATION_SPACE_REF_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

#define DEFINE_DEFAULT_CLONE_CONFIGURATION_SPACE_REF(Derived) \
  std::shared_ptr<ConfigurationSpaceRefBase> clone() const override { return std::make_shared<Derived>(*this); } 

///
/// @class ConfigurationSpaceRefBase 
/// @brief Base class of the configuration. 
///
class ConfigurationSpaceRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  ConfigurationSpaceRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~ConfigurationSpaceRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ConfigurationSpaceRefBase(const ConfigurationSpaceRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConfigurationSpaceRefBase& operator=(const ConfigurationSpaceRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConfigurationSpaceRefBase(ConfigurationSpaceRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConfigurationSpaceRefBase& operator=(ConfigurationSpaceRefBase&&) noexcept = default;

  ///
  /// @brief Clones this to a shared ptr. 
  ///
  virtual std::shared_ptr<ConfigurationSpaceRefBase> clone() const = 0;

  ///
  /// @brief Computes the reference configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] grid_info Grid info.
  /// @param[in] q_ref Reference position. Size is Robot::dimq().
  ///
  virtual void updateRef(const Robot& robot, const GridInfo& grid_info,
                         Eigen::VectorXd& q_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] grid_info Grid info.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const GridInfo& grid_info) const = 0;
};

} // namespace robotoc

#endif // ROBOTOC_CONFIGURATION_SPACE_REF_HPP_