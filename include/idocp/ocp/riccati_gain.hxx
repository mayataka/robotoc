#ifndef IDOCP_RICCATI_GAIN_HXX_
#define IDOCP_RICCATI_GAIN_HXX_

#include "idocp/ocp/riccati_gain.hpp"

#include "Eigen/LU"

#include <assert.h>


namespace idocp {

inline RiccatiGain::RiccatiGain(const Robot& robot) 
  : K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
    Ginv(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())),
    k(Eigen::VectorXd::Zero(robot.dimu())),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()) {
}


inline RiccatiGain::RiccatiGain() 
  : K(),
    Ginv(),
    k(),
    dimv_(0),
    dimu_(0) {
}


inline RiccatiGain::~RiccatiGain() {
}


inline void RiccatiGain::computeFeedbackGainAndFeedforward(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual) {
  Ginv = kkt_matrix.Quu().llt().solve(Eigen::MatrixXd::Identity(dimu_, dimu_));
  K.noalias() = - Ginv * kkt_matrix.Qxu().transpose();
  k.noalias() = - Ginv * kkt_residual.lu();
}


inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kq() const {
  return K.topLeftCorner(dimu_, dimv_);
}

inline const Eigen::Block<const Eigen::MatrixXd> RiccatiGain::Kv() const {
  return K.topRightCorner(dimu_, dimv_);
}

} // namespace idocp 

#endif // IDOCP_RICCATI_GAIN_HXX_