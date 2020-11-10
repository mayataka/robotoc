#ifndef IDOCP_RICCATI_FACTORIZATION_HXX_
#define IDOCP_RICCATI_FACTORIZATION_HXX_

#include "idocp/ocp/riccati_factorization.hxx"


namespace idocp {

inline RiccatiFactorization::RiccatiFactorization(const Robot& robot) 
  : Pqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pqv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    sq(Eigen::VectorXd::Zero(robot.dimv())),
    sv(Eigen::VectorXd::Zero(robot.dimv())),
    K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
    k(Eigen::VectorXd::Zero(robot.dimu())),
    ApBK(),
    apBk(),
    BGinvBt(),
    Pi(),
    N(),
    pi(),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()) {
  if (robot.max_point_contacts() > 0) {
    ApBK.resize(2*robot.dimv(), 2*robot.dimv());
    ApBK.setZero();
    apBk.resize(2*robot.dimv());
    apBk.setZero();
    BGinvBt.resize(robot.dimv(), robot.dimv());
    BGinvBt.setZero();
    Pi.resize(2*robot.dimv(), 2*robot.dimv());
    Pi.setZero();
    N.resize(2*robot.dimv(), 2*robot.dimv());
    N.setZero();
    pi.resize(2*robot.dimv());
    pi.setZero();
  }
}


inline RiccatiFactorization::RiccatiFactorization() 
  : Pqq(),
    Pqv(),
    Pvq(),
    Pvv(),
    sq(),
    sv(),
    K(),
    k(),
    ApBK(),
    apBk(),
    BGinvBt(),
    Pi(),
    N(),
    pi(),
    dimv_(0),
    dimu_(0) {
}


inline RiccatiFactorization::~RiccatiFactorization() {
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiFactorization::Kq() {
  return K.topLeftCorner(dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
RiccatiFactorization::Kq() const {
  return K.topLeftCorner(dimu_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> RiccatiFactorization::Kv() {
  return K.topRightCorner(dimu_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
RiccatiFactorization::Kv() const {
  return K.topRightCorner(dimu_, dimv_);
}

} // namespace idocp 

#endif // IDOCP_RICCATI_FACTORIZATION_HXX_ 