#ifndef IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HXX_
#define IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HXX_

#include "idocp/ocp/pre_impulse_split_kkt_matrix.hpp"

#include <cassert>


namespace idocp {

inline PreImpulseSplitKKTMatrix::PreImpulseSplitKKTMatrix(const Robot& robot) 
  : Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    dimv_(robot.dimv()), 
    dimf_(0) {
}


inline PreImpulseSplitKKTMatrix::PreImpulseSplitKKTMatrix() 
  : Pq_full_(),
    dimv_(0), 
    dimf_(0) {
}


inline PreImpulseSplitKKTMatrix::~PreImpulseSplitKKTMatrix() {
}


inline void PreImpulseSplitKKTMatrix::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
}


inline Eigen::Block<Eigen::MatrixXd> PreImpulseSplitKKTMatrix::Pq() {
  return Pq_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
PreImpulseSplitKKTMatrix::Pq() const {
  return Pq_full_.topLeftCorner(dimf_, dimv_);
}


inline void PreImpulseSplitKKTMatrix::setZero() {
  Pq_full_.setZero();
}


inline int PreImpulseSplitKKTMatrix::dimf() const {
  return dimf_;
}


inline bool PreImpulseSplitKKTMatrix::isApprox(
    const PreImpulseSplitKKTMatrix& other) const {
  if (!Pq().isApprox(other.Pq())) return false;
  return true;
}


inline bool PreImpulseSplitKKTMatrix::hasNaN() const {
  if (Pq_full_.hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HXX_ 