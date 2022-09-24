#ifndef ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"


namespace robotoc {

///
/// @class SwitchingConstraintData
/// @brief Data for the switching constraint.
///
class SwitchingConstraintData {
public:
  ///
  /// @brief Constructs a switching constraint.
  /// @param[in] robot Robot model. 
  ///
  SwitchingConstraintData(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  SwitchingConstraintData();

  ///
  /// @brief Default destructor. 
  ///
  ~SwitchingConstraintData() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingConstraintData(const SwitchingConstraintData&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SwitchingConstraintData& operator=(const SwitchingConstraintData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingConstraintData(SwitchingConstraintData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingConstraintData& operator=(SwitchingConstraintData&&) noexcept = default;

  ///
  /// @brief Sets the dimension of the switching constraint.
  /// @param[in] dims The dimension of the switching constraint. Must be non-negative.
  ///
  void setDimension(const int dims);

  Eigen::VectorXd q;

  Eigen::VectorXd dq;

  Eigen::VectorXd PqT_xi;

  ///
  /// @brief Jacobian of the original contact position constraint w.r.t. q. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Pq() { 
    return Pq_full_.topLeftCorner(dims_, dimv_); 
  }

  ///
  /// @brief const version of SwitchingConstraintJacobian::Pq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Pq() const {
    return Pq_full_.topLeftCorner(dims_, dimv_); 
  }

  ///
  /// @brief Returns the dimension of the switching constraint.
  /// @return Dimension of the switching constraint.
  ///
  int dims() const { return dims_; }

private:
  Eigen::MatrixXd Pq_full_;
  int dimv_, dims_;
};

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_