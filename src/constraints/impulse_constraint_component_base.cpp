#include "idocp/constraints/impulse_constraint_component_base.hpp"

#include <stdexcept>
#include <iostream>


namespace idocp {

ImpulseConstraintComponentBase::ImpulseConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rule) 
  : barrier_(barrier),
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
    try {
    if (barrier <= 0) {
      throw std::out_of_range(
          "invalid argment: barrirer must be positive");
    }
    if (fraction_to_boundary_rule <= 0) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rule must be positive");
    }
    if (fraction_to_boundary_rule >= 1) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rule must be less than 1");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


ImpulseConstraintComponentBase::ImpulseConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rule_(0) {
}

} // namespace idocp