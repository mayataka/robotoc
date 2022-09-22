#include "robotoc/core/switching_constraint_residual.hpp"


namespace robotoc {

SwitchingConstraintResidual::SwitchingConstraintResidual(
    const Robot& robot)
  : P_full_(robot.max_dimf()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dims_(0) {
}


SwitchingConstraintResidual::SwitchingConstraintResidual() 
  : P_full_(),
    dimq_(0),
    dimv_(0), 
    dims_(0) {
}

bool SwitchingConstraintResidual::isApprox(
    const SwitchingConstraintResidual& other) const {
  assert(dims() == other.dims());
  if (P().isApprox(other.P())) 
    return true;
  else 
    return false;
}


bool SwitchingConstraintResidual::hasNaN() const {
  if (P().hasNaN()) 
    return true;
  else 
    return false;
}


void SwitchingConstraintResidual::disp(std::ostream& os) const {
  os << "SwitchingConstraintResidual:" << std::endl;
  os << "  P = " << P().transpose() << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SwitchingConstraintResidual& sc_residual) {
  sc_residual.disp(os);
  return os;
}

} // namespace robotoc 