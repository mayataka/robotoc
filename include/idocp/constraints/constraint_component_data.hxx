#ifndef IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_

#include "idocp/constraints/constraint_component_data.hpp"

#include <exception>
#include <iostream>


namespace idocp {

inline ConstraintComponentData::ConstraintComponentData(const int dimc)
  : slack(Eigen::VectorXd::Zero(dimc)),
    dual(Eigen::VectorXd::Zero(dimc)),
    residual(Eigen::VectorXd::Zero(dimc)),
    duality(Eigen::VectorXd::Zero(dimc)),
    dslack(Eigen::VectorXd::Zero(dimc)),
    ddual(Eigen::VectorXd::Zero(dimc)),
    dimc_(dimc) {
  try {
    if (dimc < 0) {
      throw std::out_of_range("invalid argment: dimc must not be negative");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ConstraintComponentData::ConstraintComponentData()
  : slack(),
    dual(),
    residual(),
    duality(),
    dslack(),
    ddual(),
    dimc_(0) {
}


inline ConstraintComponentData::~ConstraintComponentData() {
}


inline int ConstraintComponentData::dimc() const {
  return dimc_;
}


inline bool ConstraintComponentData::checkDimensionalConsistency() const {
  if (slack.size() != dimc_) {
    return false;
  }
  if (dual.size() != dimc_) {
    return false;
  }
  if (residual.size() != dimc_) {
    return false;
  }
  if (duality.size() != dimc_) {
    return false;
  }
  if (dslack.size() != dimc_) {
    return false;
  }
  if (ddual.size() != dimc_) {
    return false;
  }
  return true;
}

} // namespace idocp

#endif // IDOCP_CONSTRAINT_COMPONENT_DATA_HXX_