#ifndef IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 
#define IDOCP_SPLIT_KKT_RESIDUAL_HXX_

#include "idocp/ocp/split_kkt_residual.hpp"


namespace idocp {

inline SplitKKTResidual::SplitKKTResidual(const Robot& robot) 
  : Fx(Eigen::VectorXd::Zero(2*robot.dimv())),
    lx(Eigen::VectorXd::Zero(2*robot.dimv())),
    la(Eigen::VectorXd::Zero(robot.dimv())),
    lu(Eigen::VectorXd::Zero(robot.dimu())),
    lf_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimu_(robot.dimu()),
    dimf_(0) {
}


inline SplitKKTResidual::SplitKKTResidual() 
  : Fx(),
    lx(),
    la(),
    lu(),
    lf_full_(),
    dimv_(0), 
    dimu_(0),
    dimf_(0) {
}


inline SplitKKTResidual::~SplitKKTResidual() {
}


inline void SplitKKTResidual::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fq() {
  return Fx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fq() const {
  return Fx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::Fv() {
  return Fx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::Fv() const {
  return Fx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lq() {
  return lx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lq() const {
  return lx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lv() {
  return lx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lv() const {
  return lx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTResidual::lf() {
  return lf_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTResidual::lf() const {
  return lf_full_.head(dimf_);
}


inline double SplitKKTResidual::KKTError() const {
  double err = 0;
  err += Fx.squaredNorm();
  err += lx.squaredNorm();
  err += lu.squaredNorm();
  err += la.squaredNorm();
  err += lf().squaredNorm();
  return err;
}


inline double SplitKKTResidual::constraintViolation() const {
  return Fx.template lpNorm<1>();
}


inline void SplitKKTResidual::setZero() {
  Fx.setZero();
  lx.setZero();
  la.setZero();
  lu.setZero();
  lf().setZero();
}


inline int SplitKKTResidual::dimf() const {
  return dimf_;
}


inline bool SplitKKTResidual::isDimensionConsistent() const {
  if (Fx.size() != 2*dimv_) return false;
  if (lx.size() != 2*dimv_) return false;
  if (la.size() != dimv_) return false;
  if (lu.size() != dimu_) return false;
  return true;
}


inline bool SplitKKTResidual::isApprox(const SplitKKTResidual& other) const {
  assert(isDimensionConsistent());
  assert(other.isDimensionConsistent());
  if (!Fx.isApprox(other.Fx)) return false;
  if (!lx.isApprox(other.lx)) return false;
  if (!la.isApprox(other.la)) return false;
  if (!lu.isApprox(other.lu)) return false;
  if (dimf_ > 0) {
    if (!lf().isApprox(other.lf())) return false;
  }
  return true;
}


inline bool SplitKKTResidual::hasNaN() const {
  assert(isDimensionConsistent());
  if (Fx.hasNaN()) return true;
  if (lx.hasNaN()) return true;
  if (la.hasNaN()) return true;
  if (lu.hasNaN()) return true;
  if (lf().hasNaN()) return true;
  return false;
}


inline void SplitKKTResidual::setRandom() {
  Fx.setRandom();
  lx.setRandom();
  la.setRandom();
  lu.setRandom();
  lf().setRandom();
}


inline void SplitKKTResidual::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


inline SplitKKTResidual SplitKKTResidual::Random(const Robot& robot) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom();
  return kkt_residual;
}


inline SplitKKTResidual SplitKKTResidual::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setRandom(contact_status);
  return kkt_residual;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_RESIDUAL_HXX_ 