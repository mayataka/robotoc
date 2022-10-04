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


inline bool SplitConstrainedRiccatiFactorization::isApprox(
    const SplitConstrainedRiccatiFactorization& other) const {
  if (dims() != other.dims()) return false;
  if (!DGinv().isApprox(other.DGinv())) return false;
  if (!S().isApprox(other.S())) return false;
  if (!Sinv().isApprox(other.Sinv())) return false;
  if (!SinvDGinv().isApprox(other.SinvDGinv())) return false;
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
  if (Ginv.hasNaN()) return true;
  if (DtM.hasNaN()) return true;
  if (KtDtM.hasNaN()) return true;
  return false;
}

} // namespace robotoc

#endif // ROBOTOC_SPLIT_CONSTRAINTED_RICCATI_FACTORIZATION_HXX_ 