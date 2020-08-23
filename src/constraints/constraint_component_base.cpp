#include "idocp/constraints/constraint_component_base.hpp"

#include <assert.h>


namespace idocp {

ConstraintComponentBase::ConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rate) 
  : barrier_(barrier),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
  assert(barrier > 0);
  assert(fraction_to_boundary_rate > 0);
}


ConstraintComponentBase::ConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rate_(0) {
}

} // namespace idocp