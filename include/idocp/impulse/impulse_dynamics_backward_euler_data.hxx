#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HXX_ 
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HXX_

#include "idocp/impulse/impulse_dynamics_backward_euler_data.hpp"

namespace idocp {

inline ImpulseDynamicsBackwardEulerData::ImpulseDynamicsBackwardEulerData(
    const Robot& robot) 
  : dImDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dImDddv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Minv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Qdvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    ImD(Eigen::VectorXd::Zero(robot.dimv())),
    Minv_ImD(Eigen::VectorXd::Zero(robot.dimv())),
    ldv(Eigen::VectorXd::Zero(robot.dimv())),
    Qdvf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    dimv_(robot.dimv()),
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEulerData::ImpulseDynamicsBackwardEulerData() 
  : dImDdq(),
    dImDddv(),
    Minv(),
    Qdvq(),
    ImD(),
    Minv_ImD(),
    ldv(),
    Qdvf_full_(),
    dimv_(0),
    dimf_(0) {
}


inline ImpulseDynamicsBackwardEulerData::~ImpulseDynamicsBackwardEulerData() {
}


inline void ImpulseDynamicsBackwardEulerData::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
}


inline bool ImpulseDynamicsBackwardEulerData::checkDimensions() const {
  if (dImDdq.rows() != dimv_) return false;
  if (dImDdq.cols() != dimv_) return false;
  if (dImDddv.rows() != dimv_) return false;
  if (dImDddv.cols() != dimv_) return false;
  if (Minv.rows() != dimv_) return false;
  if (Minv.cols() != dimv_) return false;
  if (Qdvq.rows() != dimv_) return false;
  if (Qdvq.cols() != dimv_) return false;
  if (ImD.size() != dimv_) return false;
  if (Minv_ImD.size() != dimv_) return false;
  if (ldv.size() != dimv_) return false;
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsBackwardEulerData::Qdvf() {
  return Qdvf_full_.topLeftCorner(dimv_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsBackwardEulerData::Qdvf() const {
  return Qdvf_full_.topLeftCorner(dimv_, dimf_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_DATA_HXX_ 