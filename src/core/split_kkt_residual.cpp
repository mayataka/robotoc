#include "robotoc/core/split_kkt_residual.hpp"

#include <random>

namespace robotoc {

SplitKKTResidual::SplitKKTResidual(const Robot& robot) 
  : Fx(Eigen::VectorXd::Zero(2*robot.dimv())),
    lx(Eigen::VectorXd::Zero(2*robot.dimv())),
    la(Eigen::VectorXd::Zero(robot.dimv())),
    ldv(Eigen::VectorXd::Zero(robot.dimv())),
    lu(Eigen::VectorXd::Zero(robot.dimu())),
    h(0.0),
    kkt_error(0.0),
    cost(0.0),
    constraint_violation(0.0),
    P_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()),
    dimf_(0),
    dims_(0) {
}


SplitKKTResidual::SplitKKTResidual() 
  : Fx(),
    lx(),
    la(),
    ldv(),
    lu(),
    h(0.0),
    kkt_error(0.0),
    cost(0.0),
    constraint_violation(0.0),
    P_full_(),
    lf_full_(),
    dimv_(0), 
    dimu_(0),
    dimf_(0),
    dims_(0) {
}


bool SplitKKTResidual::isDimensionConsistent() const {
  if (Fx.size() != 2*dimv_) return false;
  if (lx.size() != 2*dimv_) return false;
  if (la.size() != dimv_) return false;
  if (ldv.size() != dimv_) return false;
  if (lu.size() != dimu_) return false;
  return true;
}


bool SplitKKTResidual::isApprox(const SplitKKTResidual& other) const {
  assert(isDimensionConsistent());
  assert(other.isDimensionConsistent());
  if (!Fx.isApprox(other.Fx)) return false;
  if (dims_ > 0) {
    assert(dims() == other.dims());
    if (!P().isApprox(other.P())) return false;
  }
  if (!lx.isApprox(other.lx)) return false;
  if (!la.isApprox(other.la)) return false;
  if (!ldv.isApprox(other.ldv)) return false;
  if (!lu.isApprox(other.lu)) return false;
  if (dimf_ > 0) {
    assert(dimf() == other.dimf());
    if (!lf().isApprox(other.lf())) return false;
  }
  Eigen::VectorXd vec(4), other_vec(4);
  vec << h, kkt_error, cost, constraint_violation;
  other_vec << other.h, other.kkt_error, other.cost, other.constraint_violation;
  if (!vec.isApprox(other_vec)) return false;
  return true;
}


bool SplitKKTResidual::hasNaN() const {
  assert(isDimensionConsistent());
  if (Fx.hasNaN()) return true;
  if (dims() > 0) {
    if (P().hasNaN()) return true;
  }
  if (lx.hasNaN()) return true;
  if (la.hasNaN()) return true;
  if (ldv.hasNaN()) return true;
  if (lu.hasNaN()) return true;
  if (dimf() > 0) {
    if (lf().hasNaN()) return true;
  }
  Eigen::VectorXd vec(4), other_vec(4);
  vec << h, kkt_error, cost, constraint_violation;
  if (vec.hasNaN()) return true;
  return false;
}


void SplitKKTResidual::setRandom() {
  Fx.setRandom();
  P().setRandom();
  lx.setRandom();
  la.setRandom();
  ldv.setRandom();
  lu.setRandom();
  lf().setRandom();
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(4);
  h = vec.coeff(0);
  kkt_error = std::abs(vec.coeff(1));
  cost = std::abs(vec.coeff(2));
  constraint_violation = std::abs(vec.coeff(3));
}


void SplitKKTResidual::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


void SplitKKTResidual::setRandom(const ImpulseStatus& impulse_status) {
  setContactStatus(impulse_status);
  setRandom();
}


SplitKKTResidual SplitKKTResidual::Random(const Robot& robot) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom();
  return kkt_residual;
}


SplitKKTResidual SplitKKTResidual::Random(const Robot& robot, 
                                          const ContactStatus& contact_status) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom(contact_status);
  return kkt_residual;
}


SplitKKTResidual SplitKKTResidual::Random(const Robot& robot,   
                                          const ImpulseStatus& impulse_status) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom(impulse_status);
  return kkt_residual;
}


void SplitKKTResidual::disp(std::ostream& os) const {
  os << "SplitKKTResidual:" << std::endl;
  os << "  Fq = " << Fq().transpose() << std::endl;
  os << "  Fv = " << Fv().transpose() << std::endl;
  if (dims_ > 0) {
    os << "  P = " << P().transpose() << std::endl;
  }
  os << "  lq = " << lq().transpose() << std::endl;
  os << "  lv = " << lv().transpose() << std::endl;
  os << "  lu = " << lu.transpose() << std::endl;
  os << "  la = " << la.transpose() << std::endl;
  if (dimf_ > 0) {
    os << "  lf = " << lf().transpose() << std::endl;
  }
  os << "  h = " << h << std::endl;
  os << "  kkt_error = " << kkt_error << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitKKTResidual& kkt_residual) {
  kkt_residual.disp(os);
  return os;
}

} // namespace robotoc 