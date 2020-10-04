#ifndef IDOCP_RICCATI_GAIN_HXX_
#define IDOCP_RICCATI_GAIN_HXX_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

inline RiccatiGain::RiccatiGain(const Robot& robot) 
  : K_(Eigen::MatrixXd::Zero(
          robot.dimv()+2*robot.max_dimf()+robot.dim_passive(), 
          2*robot.dimv())),
    k_(Eigen::VectorXd::Zero(
          robot.dimv()+2*robot.max_dimf()+robot.dim_passive())),
    dimv_(robot.dimv()),
    dim_passive_(robot.dim_passive()),
    dimf_(robot.dimf()),
    dimc_(robot.dim_passive()+robot.dimf()) {
}


inline RiccatiGain::RiccatiGain() 
  : K_(),
    k_(),
    dimv_(0),
    dim_passive_(0),
    dimf_(0),
    dimc_(0) {
}


inline RiccatiGain::~RiccatiGain() {
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kaq() const {
  return K_.block(0, 0, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kav() const {
  return K_.block(0, dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kfq() const {
  return K_.block(dimv_, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kfv() const {
  return K_.block(dimv_, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kmuq() const {
  return K_.block(dimv_+dimf_, 0, dimc_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kmuv() const {
  return K_.block(dimv_+dimf_, dimv_, dimc_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
RiccatiGain::ka() const {
  return k_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
RiccatiGain::kf() const {
  return k_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
RiccatiGain::kmu() const {
  return k_.segment(dimv_+dimf_, dimc_);
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void RiccatiGain::computeFeedbackGain(
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
inline void RiccatiGain::computeFeedforward(
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


inline void RiccatiGain::setContactStatus(const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = dim_passive_ + contact_status.dimf();
}

} // namespace idocp 

#endif // IDOCP_RICCATI_GAIN_HXX_