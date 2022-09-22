#ifndef ROBOTOC_SPLIT_CONSTRAINTED_RICCATI_FACTORIZATION_HXX_ 
#define ROBOTOC_SPLIT_CONSTRAINTED_RICCATI_FACTORIZATION_HXX_

#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"

#include <cassert>

namespace robotoc {

inline SplitConstrainedRiccatiFactorization::
SplitConstrainedRiccatiFactorization(const Robot& robot) 
  : Ginv(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    DtM(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    KtDtM(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    DGinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    S_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Sinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    SinvDGinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    M_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    m_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mt_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mt_next_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dims_(0) { 
}


inline SplitConstrainedRiccatiFactorization::
SplitConstrainedRiccatiFactorization() 
  : Ginv(),
    DtM(),
    KtDtM(),
    DGinv_full_(),
    S_full_(),
    Sinv_full_(),
    SinvDGinv_full_(),
    M_full_(),
    m_full_(),
    mt_full_(),
    mt_next_full_(),
    dimv_(0),
    dimx_(0),
    dimu_(0),
    dims_(0) { 
}


inline SplitConstrainedRiccatiFactorization::
~SplitConstrainedRiccatiFactorization() { 
}


inline void SplitConstrainedRiccatiFactorization::setConstraintDimension(
    const int dimi) {
  dims_ = dimi;
}


inline int SplitConstrainedRiccatiFactorization::dims() const {
  return dims_;
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::DGinv() {
  return DGinv_full_.topLeftCorner(dims_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::DGinv() const {
  return DGinv_full_.topLeftCorner(dims_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::S() {
  return S_full_.topLeftCorner(dims_, dims_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::S() const {
  return S_full_.topLeftCorner(dims_, dims_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::Sinv() {
  return Sinv_full_.topLeftCorner(dims_, dims_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::Sinv() const {
  return Sinv_full_.topLeftCorner(dims_, dims_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::SinvDGinv() {
  return SinvDGinv_full_.topLeftCorner(dims_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::SinvDGinv() const {
  return SinvDGinv_full_.topLeftCorner(dims_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::M() {
  return M_full_.topLeftCorner(dims_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitConstrainedRiccatiFactorization::M() const {
  return M_full_.topLeftCorner(dims_, dimx_);
}

inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::m() {
  return m_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::m() const {
  return m_full_.head(dims_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::mt() {
  return mt_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::mt() const {
  return mt_full_.head(dims_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::mt_next() {
  return mt_next_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitConstrainedRiccatiFactorization::mt_next() const {
  return mt_next_full_.head(dims_);
}


inline bool SplitConstrainedRiccatiFactorization::isApprox(
    const SplitConstrainedRiccatiFactorization& other) const {
  if (dims() != other.dims()) return false;
  if (!DGinv().isApprox(other.DGinv())) return false;
  if (!S().isApprox(other.S())) return false;
  if (!Sinv().isApprox(other.Sinv())) return false;
  if (!SinvDGinv().isApprox(other.SinvDGinv())) return false;
  if (!M().isApprox(other.M())) return false;
  if (!m().isApprox(other.m())) return false;
  if (!mt().isApprox(other.mt())) return false;
  if (!mt_next().isApprox(other.mt_next())) return false;
  if (!Ginv.isApprox(other.Ginv)) return false;
  if (!DtM.isApprox(other.DtM)) return false;
  if (!KtDtM.isApprox(other.KtDtM)) return false;
  return true;
}


inline bool SplitConstrainedRiccatiFactorization::hasNaN() const {
  if (DGinv().hasNaN()) return true;
  if (S().hasNaN()) return true;
  if (Sinv().hasNaN()) return true;
  if (SinvDGinv().hasNaN()) return true;
  if (M().hasNaN()) return true;
  if (m().hasNaN()) return true;
  if (mt().hasNaN()) return true;
  if (mt_next().hasNaN()) return true;
  if (Ginv.hasNaN()) return true;
  if (DtM.hasNaN()) return true;
  if (KtDtM.hasNaN()) return true;
  return false;
}

} // namespace robotoc

#endif // ROBOTOC_SPLIT_CONSTRAINTED_RICCATI_FACTORIZATION_HXX_ 