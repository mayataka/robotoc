#ifndef ROBOTOC_SPLIT_KKT_MATRIX_HXX_
#define ROBOTOC_SPLIT_KKT_MATRIX_HXX_

#include "robotoc/ocp/split_kkt_matrix.hpp"

#include <cassert>


namespace robotoc {

inline SplitKKTMatrix::SplitKKTMatrix(const Robot& robot) 
  : Fxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Fvu(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimu())),
    Qxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Qaa(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qxu(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimu())),
    Quu(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    Fqq_prev(),
    fx(Eigen::VectorXd::Zero(2*robot.dimv())),
    Qtt(0),
    Qtt_prev(0),
    hx(Eigen::VectorXd::Zero(2*robot.dimv())),
    hu(Eigen::VectorXd::Zero(robot.dimu())),
    Qff_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Qqf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()), 
    dimf_(0) {
  if (robot.hasFloatingBase()) {
    Fqq_prev.resize(robot.dimv(), robot.dimv());
    Fqq_prev.setZero();
  }
}


inline SplitKKTMatrix::SplitKKTMatrix() 
  : Fxx(),
    Fvu(),
    Qxx(),
    Qaa(),
    Qxu(),
    Quu(),
    Fqq_prev(),
    fx(),
    Qtt(0),
    Qtt_prev(0),
    hx(),
    hu(),
    Qff_full_(),
    Qqf_full_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimu_(0), 
    dimf_(0) {
}


inline SplitKKTMatrix::~SplitKKTMatrix() {
}


inline void SplitKKTMatrix::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqq() {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqq() const {
  return Fxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fqv() {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fqv() const {
  return Fxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvq() {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvq() const {
  return Fxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Fvv() {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Fvv() const {
  return Fxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqq() {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqq() const {
  return Qxx.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqv() {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqv() const {
  return Qxx.topRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvq() {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvq() const {
  return Qxx.bottomLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvv() {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvv() const {
  return Qxx.bottomRightCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqu() {
  return Qxu.topLeftCorner(dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqu() const {
  return Qxu.topLeftCorner(dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qvu() {
  return Qxu.bottomLeftCorner(dimv_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qvu() const {
  return Qxu.bottomLeftCorner(dimv_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qff() {
  return Qff_full_.topLeftCorner(dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qff() const {
  return Qff_full_.topLeftCorner(dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitKKTMatrix::Qqf() {
  return Qqf_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitKKTMatrix::Qqf() const {
  return Qqf_full_.topLeftCorner(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::fq() {
  return fx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTMatrix::fq() const {
  return fx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::fv() {
  return fx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTMatrix::fv() const {
  return fx.tail(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::hq() {
  return hx.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTMatrix::hq() const {
  return hx.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitKKTMatrix::hv() {
  return hx.tail(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitKKTMatrix::hv() const {
  return hx.tail(dimv_);
}


inline void SplitKKTMatrix::setZero() {
  Fxx.setZero();
  Fvu.setZero();
  Qxx.setZero();
  Qaa.setZero();
  Qxu.setZero();
  Quu.setZero();
  Qff().setZero();
  Qqf().setZero();
  Fqq_prev.setZero();
  fx.setZero();
  Qtt = 0;
  Qtt_prev = 0;
  hx.setZero();
  hu.setZero();
}


inline int SplitKKTMatrix::dimf() const {
  return dimf_;
}


inline bool SplitKKTMatrix::isDimensionConsistent() const {
  if (Fxx.rows() != 2*dimv_) return false;
  if (Fxx.cols() != 2*dimv_) return false;
  if (Fvu.rows() != dimv_) return false;
  if (Fvu.cols() != dimu_) return false;
  if (Qxx.rows() != 2*dimv_) return false;
  if (Qxx.cols() != 2*dimv_) return false;
  if (Qaa.rows() != dimv_) return false;
  if (Qaa.cols() != dimv_) return false;
  if (Qxu.rows() != 2*dimv_) return false;
  if (Qxu.cols() != dimu_) return false;
  if (Quu.rows() != dimu_) return false;
  if (Quu.cols() != dimu_) return false;
  if (has_floating_base_) {
    if (Fqq_prev.rows() != dimv_) return false;
    if (Fqq_prev.cols() != dimv_) return false;
  }
  if (fx.size() != 2*dimv_) return false;
  if (hx.size() != 2*dimv_) return false;
  if (hu.size() != dimu_) return false;
  return true;
}


inline bool SplitKKTMatrix::isApprox(const SplitKKTMatrix& other) const {
  if (!Fxx.isApprox(other.Fxx)) return false;
  if (!Fvu.isApprox(other.Fvu)) return false;
  if (!Qxx.isApprox(other.Qxx)) return false;
  if (!Qaa.isApprox(other.Qaa)) return false;
  if (!Qxu.isApprox(other.Qxu)) return false;
  if (!Quu.isApprox(other.Quu)) return false;
  if (!Qff().isApprox(other.Qff())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!Fqq_prev.isApprox(other.Fqq_prev)) return false;
  if (!fx.isApprox(other.fx)) return false;
  Eigen::VectorXd vec(2), other_vec(2);
  vec << Qtt, Qtt_prev;
  other_vec << other.Qtt, other.Qtt_prev;
  if (!vec.isApprox(other_vec)) return false;
  if (!hx.isApprox(other.hx)) return false;
  if (!hu.isApprox(other.hu)) return false;
  return true;
}


inline bool SplitKKTMatrix::hasNaN() const {
  if (Fxx.hasNaN()) return true;
  if (Fvu.hasNaN()) return true;
  if (Qxx.hasNaN()) return true;
  if (Qaa.hasNaN()) return true;
  if (Qxu.hasNaN()) return true;
  if (Quu.hasNaN()) return true;
  if (Qff().hasNaN()) return true;
  if (Qqf().hasNaN()) return true;
  if (Fqq_prev.hasNaN()) return true;
  if (fx.hasNaN()) return true;
  Eigen::VectorXd vec(2);
  vec << Qtt, Qtt_prev;
  if (vec.hasNaN()) return true;
  if (hx.hasNaN()) return true;
  if (hu.hasNaN()) return true;
  return false;
}


inline void SplitKKTMatrix::setRandom() {
  Fxx.setRandom();
  Fvu.setRandom();
  const Eigen::MatrixXd Qxxuu_seed = Eigen::MatrixXd::Random(dimx_+dimu_, dimx_+dimu_);
  const Eigen::MatrixXd Qxxuu = Qxxuu_seed * Qxxuu_seed.transpose();
  Qxx = Qxxuu.topLeftCorner(dimx_, dimx_);
  Qxu = Qxxuu.topRightCorner(dimx_, dimu_);
  Quu = Qxxuu.bottomRightCorner(dimu_, dimu_);
  const Eigen::MatrixXd Qaaff_seed = Eigen::MatrixXd::Random(dimv_+dimf_, dimv_+dimf_);
  const Eigen::MatrixXd Qaaff = Qaaff_seed * Qaaff_seed.transpose();
  Qaa = Qaaff.topLeftCorner(dimv_, dimv_);
  Qff() = Qaaff.bottomRightCorner(dimf_, dimf_);
  Qqf().setRandom();
  Fqq_prev.setRandom();
  fx.setRandom();
  Qtt = Eigen::VectorXd::Random(1)[0];
  Qtt_prev = Eigen::VectorXd::Random(1)[0];
  hx.setRandom();
  hu.setRandom();
}


inline void SplitKKTMatrix::setRandom(const ContactStatus& contact_status) {
  setContactStatus(contact_status);
  setRandom();
}


inline SplitKKTMatrix SplitKKTMatrix::Random(const Robot& robot) {
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom();
  return kkt_matrix;
}


inline SplitKKTMatrix SplitKKTMatrix::Random(
    const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setRandom(contact_status);
  return kkt_matrix;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_KKT_MATRIX_HXX_ 