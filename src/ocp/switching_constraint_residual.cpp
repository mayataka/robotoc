#include "robotoc/ocp/switching_constraint_residual.hpp"


namespace robotoc {

void SwitchingConstraintResidual::disp(std::ostream& os) const {
  os << "switching constraint residual:" << std::endl;
  os << "  P = " << P().transpose() << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SwitchingConstraintResidual& sc_residual) {
  sc_residual.disp(os);
  return os;
}

} // namespace robotoc 