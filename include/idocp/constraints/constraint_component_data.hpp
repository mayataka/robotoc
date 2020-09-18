#ifndef IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_

#include "Eigen/Core"


namespace idocp {

///
/// @class ConstraintComponentData
/// @brief Data used in constraint components. Composed by slack, 
/// dual (Lagrange multiplier), primal residual, duality between the slack and 
/// dual, and directions of slack and dual.
///
class ConstraintComponentData {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] dimc Dimension of the constraint component. 
  ///
  ConstraintComponentData(const int dimc);

  ///
  /// @brief Default constructor. 
  ///
  ConstraintComponentData();

  ///
  /// @brief Destructor. 
  ///
  ~ConstraintComponentData();

  ///
  /// @brief Default copy constructor. 
  ///
  ConstraintComponentData(const ConstraintComponentData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConstraintComponentData& operator=(const ConstraintComponentData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConstraintComponentData(ConstraintComponentData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConstraintComponentData& operator=(ConstraintComponentData&&) noexcept 
      = default;

  ///
  /// @brief Slack variable of the constraint. Size is 
  /// ConstraintComponentData::dimc(). All elements must be positive.
  ///
  Eigen::VectorXd slack;

  ///
  /// @brief Dual variable (Lagrange multiplier) of the constraint. Size is 
  /// ConstraintComponentData::dimc(). All elements must be positive.
  ///
  Eigen::VectorXd dual;

  ///
  /// @brief Primal residual of the constraint. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd residual;

  ///
  /// @brief Residual of the duality between slakc and dual. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd duality;

  ///
  /// @brief Newton direction of the slack. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd dslack;

  ///
  /// @brief Newton direction of the dual. Size is 
  /// ConstraintComponentData::dimc(). 
  ///
  Eigen::VectorXd ddual;

  ///
  /// @brief Dimension of the constraint. 
  /// @return Dimension of the constraint. 
  ///
  int dimc() const;

  ///
  /// @brief Check whether dimensions of slack, dual, residual, duality, 
  /// dslack, ddual are ConstraintComponentData::dimc(). 
  /// @return Dimension of the constraint. 
  ///
  bool checkDimensionnalConsistency() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimc_;

};

} // namespace idocp

#include "idocp/constraints/constraint_component_data.hxx"

#endif // IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_