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
    M_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    m_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mt_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    mt_next_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dims_(0) {
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
    M_full_(),
    m_full_(),
    mt_full_(),
    mt_next_full_(),
    dimv_(0),
    dimx_(0),
    dims_(0) {
}


inline SplitRiccatiFactorization::~SplitRiccatiFactorization() {
}


inline void SplitRiccatiFactorization::setConstraintDimension(const int dims) {
  assert(dims >= 0);
  assert(dims <= m_full_.size());
  dims_ = dims;
}


inline int SplitRiccatiFactorization::dims() const {
  return dims_;
}


inline Eigen::Block<Eigen::MatrixXd> SplitRiccatiFactorization::M() {
  return M_full_.topLeftCorner(dims_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> SplitRiccatiFactorization::M() const {
  return M_full_.topLeftCorner(dims_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitRiccatiFactorization::m() {
  return m_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitRiccatiFactorization::m() const {
  return m_full_.head(dims_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitRiccatiFactorization::mt() {
  return mt_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitRiccatiFactorization::mt() const {
  return mt_full_.head(dims_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitRiccatiFactorization::mt_next() {
  return mt_next_full_.head(dims_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> SplitRiccatiFactorization::mt_next() const {
  return mt_next_full_.head(dims_);
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
  M_full_.setZero();
  m_full_.setZero();
  mt_full_.setZero();
  mt_next_full_.setZero();
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
  M_full_.setRandom();
  m_full_.setRandom();
  mt_full_.setRandom();
  mt_next_full_.setRandom();
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
  if (dims() > 0) {
    if (!M().isApprox(other.M())) return false;
    if (!m().isApprox(other.m())) return false;
    if (!mt().isApprox(other.mt())) return false;
    if (!mt_next().isApprox(other.mt_next())) return false;
  }
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
  if (dims() > 0) {
    if (M().hasNaN()) return true;
    if (m().hasNaN()) return true;
    if (mt().hasNaN()) return true;
    if (mt_next().hasNaN()) return true;
  }
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