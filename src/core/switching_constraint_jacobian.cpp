#include "robotoc/core/switching_constraint_jacobian.hpp"


namespace robotoc {

SwitchingConstraintJacobian::SwitchingConstraintJacobian(const Robot& robot)
  : Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phix_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    Phia_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phiu_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    Phit_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimi_(0) {
}


SwitchingConstraintJacobian::SwitchingConstraintJacobian()
  : Pq_full_(),
    Phix_full_(),
    Phia_full_(),
    Phiu_full_(),
    Phit_full_(),
    has_floating_base_(false),
    dimv_(0),
    dimx_(0),
    dimu_(0),
    dimi_(0) {
}


bool SwitchingConstraintJacobian::isApprox(
    const SwitchingConstraintJacobian& other) const {
  assert(dimi() == other.dimi());
  if (!Pq().isApprox(other.Pq())) return false;
  if (!Phix().isApprox(other.Phix())) return false;
  if (!Phia().isApprox(other.Phia())) return false;
  if (!Phiu().isApprox(other.Phiu())) return false;
  if (!Phit().isApprox(other.Phit())) return false;
  return true;
}


bool SwitchingConstraintJacobian::hasNaN() const {
  if (Pq().hasNaN()) return true;
  if (Phix().hasNaN()) return true;
  if (Phia().hasNaN()) return true;
  if (Phiu().hasNaN()) return true;
  if (Phit().hasNaN()) return true;
  return false;
}


void SwitchingConstraintJacobian::disp(std::ostream& os) const {
  os << "SwitchingConstraintJacobian:" << std::endl;
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