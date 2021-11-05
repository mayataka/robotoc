#ifndef ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_
#define ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"


namespace robotoc {

///
/// @class SplitRiccatiFactorization
/// @brief Riccati factorization matrix and vector for a time stage.
///
class SplitRiccatiFactorization {
public:
  ///
  /// @brief Constructs Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. 
  ///
  SplitRiccatiFactorization(const Robot& robot)
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

  ///
  /// @brief Default constructor. 
  ///
  SplitRiccatiFactorization()
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

  ///
  /// @brief Destructor. 
  ///
  ~SplitRiccatiFactorization() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  SplitRiccatiFactorization(const SplitRiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitRiccatiFactorization& operator=(
      const SplitRiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitRiccatiFactorization(
      SplitRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitRiccatiFactorization& operator=(
      SplitRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Riccati factorization matrix. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::MatrixXd P;

  ///
  /// @brief Riccati factorization vector. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd s;

  Eigen::Block<Eigen::MatrixXd> Pqq() {
    return P.topLeftCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pqq() const {
    return P.topLeftCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pqv() {
    return P.topRightCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pqv() const {
    return P.topRightCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pvq() {
    return P.bottomLeftCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pvq() const {
    return P.bottomLeftCorner(dimv_, dimv_); 
  }

  Eigen::Block<Eigen::MatrixXd> Pvv() {
    return P.bottomRightCorner(dimv_, dimv_); 
  }

  const Eigen::Block<const Eigen::MatrixXd> Pvv() const {
    return P.bottomRightCorner(dimv_, dimv_); 
  }

  Eigen::VectorBlock<Eigen::VectorXd> sq() {
    return s.head(dimv_);
  }

  const Eigen::VectorBlock<const Eigen::VectorXd> sq() const {
    return s.head(dimv_);
  }

  Eigen::VectorBlock<Eigen::VectorXd> sv() {
    return s.tail(dimv_);
  }

  const Eigen::VectorBlock<const Eigen::VectorXd> sv() const {
    return s.tail(dimv_);
  }

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd psi_x;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// Robot::dimu().
  ///
  Eigen::VectorXd psi_u;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd Psi;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd phi_x;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// Robot::dimu().
  ///
  Eigen::VectorXd phi_u;

  ///
  /// @brief Riccati factorization vector w.r.t. the switching time. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorXd Phi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double xi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double chi;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double rho;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double eta;

  ///
  /// @brief Riccati factorization w.r.t. the switching time. 
  ///
  double iota;

  void setZero() {
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

  void setRandom() {
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

  ///
  /// @brief Checks the equivalence of two SplitRiccatiFactorization.
  /// @param[in] other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitRiccatiFactorization& other) const {
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

  ///
  /// @brief Checks this object has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const {
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

  static SplitRiccatiFactorization Random(const Robot& robot) {
    auto riccati = SplitRiccatiFactorization(robot);
    riccati.setRandom();
    return riccati;
  }

private:
  int dimv_, dimx_;

};

} // namespace robotoc 

#endif // ROBOTOC_SPLIT_RICCATI_FACTORIZATION_HPP_ 