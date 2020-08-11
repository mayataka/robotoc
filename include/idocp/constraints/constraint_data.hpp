#ifndef IDOCP_CONSTRAINT_DATA_HPP_
#define IDOCP_CONSTRAINT_DATA_HPP_

#include <vector>

#include "Eigen/Core"


namespace idocp {

class ConstraintData {
public:
  ConstraintData(const int dimc)
    : slack(Eigen::VectorXd::Zero(dimc)),
      dual(Eigen::VectorXd::Zero(dimc)),
      residual(Eigen::VectorXd::Zero(dimc)),
      duality(Eigen::VectorXd::Zero(dimc)),
      dslack(Eigen::VectorXd::Zero(dimc)),
      ddual(Eigen::VectorXd::Zero(dimc)) {
  }

  ConstraintData()
    : slack(),
      dual(),
      residual(),
      duality(),
      dslack(),
      ddual() {
  }

  ~ConstraintData() {
  }

  // Use default copy constructor.
  ConstraintData(const ConstraintData&) = default;

  // Use default copy coperator.
  ConstraintData& operator=(const ConstraintData&) = default;

  // Use default move constructor.
  ConstraintData(ConstraintData&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintData& operator=(ConstraintData&&) noexcept = default;

  Eigen::VectorXd slack, dual, residual, duality, dslack, ddual;

};

} // namespace idocp


#endif // IDOCP_CONSTRAINT_DATA_HPP_