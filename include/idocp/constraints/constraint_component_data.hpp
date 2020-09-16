#ifndef IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_

#include "Eigen/Core"


namespace idocp {

class ConstraintComponentData {
public:
  ConstraintComponentData(const int dimc);

  ConstraintComponentData();

  ~ConstraintComponentData();

  // Use default copy constructor.
  ConstraintComponentData(const ConstraintComponentData&) = default;

  // Use default copy coperator.
  ConstraintComponentData& operator=(const ConstraintComponentData&) = default;

  // Use default move constructor.
  ConstraintComponentData(ConstraintComponentData&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintComponentData& operator=(ConstraintComponentData&&) noexcept 
      = default;

  Eigen::VectorXd slack, dual, residual, duality, dslack, ddual;

  int dimc() const;

  bool checkDimensionsOfVectors() const;

private:
  int dimc_;

};

} // namespace idocp

#include "idocp/constraints/constraint_component_data.hxx"

#endif // IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_