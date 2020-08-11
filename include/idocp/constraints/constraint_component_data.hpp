#ifndef IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_

#include <vector>

#include "Eigen/Core"


namespace idocp {

class ConstraintComponentData {
public:
  ConstraintComponentData(const int dimc)
    : slack(Eigen::VectorXd::Zero(dimc)),
      dual(Eigen::VectorXd::Zero(dimc)),
      residual(Eigen::VectorXd::Zero(dimc)),
      duality(Eigen::VectorXd::Zero(dimc)),
      dslack(Eigen::VectorXd::Zero(dimc)),
      ddual(Eigen::VectorXd::Zero(dimc)) {
  }

  ConstraintComponentData()
    : slack(),
      dual(),
      residual(),
      duality(),
      dslack(),
      ddual() {
  }

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

};

} // namespace idocp


#endif // IDOCP_CONSTRAINT_COMPONENT_DATA_HPP_