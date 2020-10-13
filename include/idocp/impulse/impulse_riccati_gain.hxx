#ifndef IDOCP_IMPULSE_RICCATI_GAIN_HXX_
#define IDOCP_IMPULSE_RICCATI_GAIN_HXX_

#include "idocp/impulse/impulse_riccati_gain.hpp"

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

inline ImpulseRiccatiGain::ImpulseRiccatiGain(const Robot& robot) 
  : K_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    k_(Eigen::VectorXd::Zero(3*robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0),
    dimc_(0) {
}


inline ImpulseRiccatiGain::ImpulseRiccatiGain() 
  : K_(),
    k_(),
    dimv_(0),
    dimf_(0),
    dimc_(0) {
}


inline ImpulseRiccatiGain::~ImpulseRiccatiGain() {
}


inline void ImpulseRiccatiGain::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = contact_status.dimf();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseRiccatiGain::Kfq() const {
  return K_.block(0, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseRiccatiGain::Kfv() const {
  return K_.block(0, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseRiccatiGain::Kmuq() const {
  return K_.block(dimf_, 0, dimc_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseRiccatiGain::Kmuv() const {
  return K_.block(dimf_, dimv_, dimc_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseRiccatiGain::kf() const {
  return k_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseRiccatiGain::kmu() const {
  return k_.segment(dimf_, dimc_);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void ImpulseRiccatiGain::computeFeedbackGain(
    const Eigen::MatrixBase<MatrixType1>& Ginv, 
    const Eigen::MatrixBase<MatrixType2>& Qafqv, 
    const Eigen::MatrixBase<MatrixType3>& Cqv) {
  const int dimaf = dimv_ + dimf_;
  assert(Ginv.rows() == dimaf+dimc_);
  assert(Ginv.cols() == dimaf+dimc_);
  assert(Qafqv.rows() == dimaf);
  assert(Qafqv.cols() == 2*dimv_);
  assert(Cqv.rows() == dimc_);
  assert(Cqv.cols() == 2*dimv_);
  K_.topRows(dimaf+dimc_).noalias() = - Ginv.leftCols(dimaf) * Qafqv;
  K_.topRows(dimaf+dimc_).noalias() -= Ginv.rightCols(dimc_) * Cqv;
}


template <typename MatrixType, typename VectorType1, typename VectorType2>
inline void ImpulseRiccatiGain::computeFeedforward(
    const Eigen::MatrixBase<MatrixType>& Ginv, 
    const Eigen::MatrixBase<VectorType1>& laf, 
    const Eigen::MatrixBase<VectorType2>& C) {
  const int dimaf = dimv_ + dimf_;
  assert(Ginv.rows() == dimaf+dimc_);
  assert(Ginv.cols() == dimaf+dimc_);
  assert(laf.size() == dimaf);
  assert(C.size() == dimc_);
  k_.head(dimaf+dimc_).noalias() = - Ginv.leftCols(dimaf) * laf;
  k_.head(dimaf+dimc_).noalias() -= Ginv.rightCols(dimc_) * C;
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_RICCATI_GAIN_HXX_ 