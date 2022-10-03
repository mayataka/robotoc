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
  Eigen::VectorXd vec(1), other_vec(1);
  vec << h;
  other_vec << other.h;
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
  Eigen::VectorXd vec(1);
  vec << h;
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
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(1);
  h = vec.coeff(0);
}


void SplitKKTResidual::setRandom(const ContactStatus& contact_status) {
  setContactDimension(contact_status.dimf());
  setRandom();
}


void SplitKKTResidual::setRandom(const ImpactStatus& impact_status) {
  setContactDimension(impact_status.dimf());
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
                                          const ImpactStatus& impact_status) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom(impact_status);
  return kkt_residual;
}


void SplitKKTResidual::disp(std::ostream& os) const {
  os << "SplitKKTResidual:" << "\n";
  os << "  Fq = " << Fq().transpose() << "\n";
  os << "  Fv = " << Fv().transpose() << "\n";
  if (dims() > 0) {
    os << "  P = " << P().transpose() << "\n";
  }
  os << "  lq = " << lq().transpose() << "\n";
  os << "  lv = " << lv().transpose() << "\n";
  os << "  lu = " << lu.transpose() << "\n";
  os << "  la = " << la.transpose() << "\n";
  os << "  ldv = " << ldv.transpose() << "\n";
  if (dimf() > 0) {
    os << "  lf = " << lf().transpose() << "\n";
  }
  os << "  h = " << h << "\n";
}


std::ostream& operator<<(std::ostream& os, 
                         const SplitKKTResidual& kkt_residual) {
  kkt_residual.disp(os);
  return os;
}

} // namespace robotoc 