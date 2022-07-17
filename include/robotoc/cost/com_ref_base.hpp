#ifndef ROBOTOC_COM_REF_BASE_HPP_
#define ROBOTOC_COM_REF_BASE_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

#define DEFINE_DEFAULT_CLONE_COM_REF(Derived) \
  std::shared_ptr<CoMRefBase> clone() const override { return std::make_shared<Derived>(*this); } 

///
/// @class CoMRefBase
/// @brief Base class of reference position of the center of mass (CoM). 
///
class CoMRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  CoMRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~CoMRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  CoMRefBase(const CoMRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CoMRefBase& operator=( const CoMRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CoMRefBase(CoMRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CoMRefBase& operator=(CoMRefBase&&) noexcept = default;

  ///
  /// @brief Clones this to a shared ptr. 
  ///
  virtual std::shared_ptr<CoMRefBase> clone() const = 0;

  ///
  /// @brief Computes the reference CoM position. 
  /// @param[in] grid_info Grid info.
  /// @param[in] com_ref Reference CoM position Size is 3.
  ///
  virtual void updateRef(const GridInfo& grid_info, 
                         Eigen::VectorXd& com_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not for the given grid info. 
  /// @param[in] grid_info Grid info.
  /// @return true if the cost is active for the given grid_info. false if not.
  ///
  virtual bool isActive(const GridInfo& grid_info) const = 0;
};

} // namespace robotoc

#endif // ROBOTOC_COM_REF_BASE_HPP_