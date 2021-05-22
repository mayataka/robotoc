#ifndef IDOCP_SPLIT_KKT_MATRIX_HXX_
#define IDOCP_SPLIT_KKT_MATRIX_HXX_

#include "idocp/ocp/split_kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline SplitKKTMatrix::SplitKKTMatrix(const Robot& robot) 
  : Fxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Fvu(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimu())),
    Qxx(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Qaa(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qxu(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimu())),
    Qxu_passive(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dim_passive())),
    Quu(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    Quu_passive_topRight(Eigen::MatrixXd::Zero(robot.dim_passive(), robot.dimu())),
    Fqq_prev(),
    Fqq_inv(),
    Fqq_prev_inv(),
    Qff_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Qqf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimu_(robot.dimu()), 
    dim_passive_(robot.dim_passive()),
    dimf_(0) {
  if (robot.hasFloatingBase()) {
    Fqq_prev.resize(robot.dimv(), robot.dimv());
    Fqq_prev.setZero();
    Fqq_inv.resize(6, 6);
    Fqq_inv.setZero();
    Fqq_prev_inv.resize(6, 6);
    Fqq_prev_inv.setZero();
  }
}


inline SplitKKTMatrix::SplitKKTMatrix() 
  : Fxx(),
    Fvu(),
    Qxx(),
    Qaa(),
    Qxu(),
    Quu(),
    Quu_passive_topRight(),
    Fqq_prev(),
    Fqq_inv(),
    Fqq_prev_inv(),
    Qff_full_(),
    Qqf_full_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimu_(0), 
    dim_passive_(0),
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


inline void SplitKKTMatrix::setZero() {
  Fxx.setZero();
  Fvu.setZero();
  Qxx.setZero();
  Qaa.setZero();
  Qxu.setZero();
  Qxu_passive.setZero();
  Quu.setZero();
  Quu_passive_topRight.setZero();
  Qff().setZero();
  Qqf().setZero();
  Fqq_prev.setZero();
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
  if (Qxu_passive.rows() != 2*dimv_) return false;
  if (Qxu_passive.cols() != dim_passive_) return false;
  if (Quu.rows() != dimu_) return false;
  if (Quu.cols() != dimu_) return false;
  if (Quu_passive_topRight.rows() != dim_passive_) return false;
  if (Quu_passive_topRight.cols() != dimu_) return false;
  if (has_floating_base_) {
    if (Fqq_prev.rows() != dimv_) return false;
    if (Fqq_prev.cols() != dimv_) return false;
  }
  return true;
}


inline bool SplitKKTMatrix::isApprox(const SplitKKTMatrix& other) const {
  if (!Fxx.isApprox(other.Fxx)) return false;
  if (!Fvu.isApprox(other.Fvu)) return false;
  if (!Qxx.isApprox(other.Qxx)) return false;
  if (!Qaa.isApprox(other.Qaa)) return false;
  if (!Qxu.isApprox(other.Qxu)) return false;
  if (has_floating_base_) {
    if (!Qxu_passive.isApprox(other.Qxu_passive)) return false;
  }
  if (!Quu.isApprox(other.Quu)) return false;
  if (has_floating_base_) {
    if (!Quu_passive_topRight.isApprox(other.Quu_passive_topRight)) 
        return false;
  }
  if (!Qff().isApprox(other.Qff())) return false;
  if (!Qqf().isApprox(other.Qqf())) return false;
  if (!Fqq_prev.isApprox(other.Fqq_prev)) return false;
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
  return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_KKT_MATRIX_HXX_ 