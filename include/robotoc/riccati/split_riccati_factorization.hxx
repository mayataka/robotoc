#ifndef ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HXX_
#define ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HXX_

#include "robotoc/riccati/split_riccati_factorization.hpp"


namespace robotoc {

inline SplitRiccatiFactorization::SplitRiccatiFactorization(const Robot& robot)
  : P(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    s(Eigen::VectorXd::Zero(2*robot.dimv())),
    psi_x(Eigen::VectorXd::Zero(2*robot.dimv())),
    psi_u(Eigen::VectorXd::Zero(robot.dimu())),
    Psi(Eigen::VectorXd::Zero(2*robot.dimv())),
    phi_x(Eigen::VectorXd::Zero(2*robot.dimv())),
    phi_u(Eigen::VectorXd::Zero(robot.dimu())),
    Phi(Eigen::VectorXd::Zero(2*robot.dimv())),
    xi(0.0),
    chi(0.0),
    rho(0.0),
    eta(0.0),
    iota(0.0),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()) {
}


inline SplitRiccatiFactorization::SplitRiccatiFactorization()
  : P(),
    s(),
    psi_x(),
    psi_u(),
    Psi(),
    phi_x(),
    phi_u(),
    Phi(),
    xi(0.0),
    chi(0.0),
    rho(0.0),
    eta(0.0),
    iota(0.0),
    dimv_(0.0),
    dimx_(0.0) {
}


inline SplitRiccatiFactorization::~SplitRiccatiFactorization() {
}


inline void SplitRiccatiFactorization::setZero() {
  P.setZero();
  s.setZero();
  psi_x.setZero();
  psi_u.setZero();
  Psi.setZero();
  phi_x.setZero();
  phi_u.setZero();
  Phi.setZero();
  xi = 0;
  chi = 0;
  rho = 0;
  eta = 0;
  iota = 0;
}


inline void SplitRiccatiFactorization::setRandom() {
  P.setRandom();
  s.setRandom();
  psi_x.setRandom();
  psi_u.setRandom();
  Psi.setRandom();
  phi_x.setRandom();
  phi_u.setRandom();
  Phi.setRandom();
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(5);
  xi = vec[0];
  chi = vec[1];
  rho = vec[2];
  eta = vec[3];
  iota = vec[4];
}


inline bool SplitRiccatiFactorization::isApprox(
  const SplitRiccatiFactorization& other) const {
  if (!P.isApprox(other.P)) return false;
  if (!s.isApprox(other.s)) return false;
  if (!psi_x.isApprox(other.psi_x)) return false;
  if (!psi_u.isApprox(other.psi_u)) return false;
  if (!Psi.isApprox(other.Psi)) return false;
  if (!phi_x.isApprox(other.phi_x)) return false;
  if (!phi_u.isApprox(other.phi_u)) return false;
  if (!Phi.isApprox(other.Phi)) return false;
  Eigen::VectorXd vec(5), other_vec(5);
  vec << xi, chi, rho, eta, iota;
  other_vec << other.xi, other.chi, other.rho, other.eta, other.iota;
  if (!vec.isApprox(other_vec)) return false;
  return true;
}


inline bool SplitRiccatiFactorization::hasNaN() const {
  if (P.hasNaN()) return true;
  if (s.hasNaN()) return true;
  if (psi_x.hasNaN()) return true;
  if (psi_u.hasNaN()) return true;
  if (Psi.hasNaN()) return true;
  if (phi_x.hasNaN()) return true;
  if (phi_u.hasNaN()) return true;
  if (Phi.hasNaN()) return true;
  Eigen::VectorXd vec(5);
  vec << xi, chi, rho, eta, iota;
  if (vec.hasNaN()) return true;
  return false;
}


inline SplitRiccatiFactorization SplitRiccatiFactorization::Random(
    const Robot& robot) {
  auto riccati = SplitRiccatiFactorization(robot);
  riccati.setRandom();
  return riccati;
}

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HXX_ 