#include "robotoc/ocp/switching_constraint_jacobian.hpp"


namespace robotoc {

void SwitchingConstraintJacobian::disp(std::ostream& os) const {
  os << "switching constraint Jacobian:" << std::endl;
  os << "  Pq = " << Pq() << std::endl;
  os << "  Phix = " << Phix() << std::endl;
  os << "  Phia = " << Phia() << std::endl;
  os << "  Phiu = " << Phiu() << std::endl;
  os << "  Phit = " << Phit().transpose() << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SwitchingConstraintJacobian& sc_jacobian) {
  sc_jacobian.disp(os);
  return os;
}

} // namespace robotoc 