#include "idocp/constraints/constraint_component_base.hpp"

#include <stdexcept>
#include <iostream>


namespace idocp {

ConstraintComponentBase::ConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rate) 
  : barrier_(barrier),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
  try {
    if (barrier <= 0) {
      throw std::out_of_range(
          "Invalid argment: barrirer must be positive!");
    }
    if (fraction_to_boundary_rate <= 0) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rate must be positive!");
    }
    if (fraction_to_boundary_rate >= 1) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rate must be less than 1!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


ConstraintComponentBase::ConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rate_(0) {
}

} // namespace idocp